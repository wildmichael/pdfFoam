/*---------------------------------------------------------------------------*\
                pdfFoam: General Purpose PDF Solution Algorithm
                   for Reactive Flow Simulations in OpenFOAM

 Copyright (C) 2012 Michael Wild, Heng Xiao, Patrick Jenny,
                    Institute of Fluid Dynamics, ETH Zurich
-------------------------------------------------------------------------------
License
    This file is part of pdfFoam.

    pdfFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) version 3 of the same License.

    pdfFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with pdfFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mcInletOutletBoundary.H"

#include "addToRunTimeSelectionTable.H"
#include "mcParticle.H"
#include "mcParticleCloud.H"
#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "Random.H"
#include "transformField.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "DynamicList.H"
#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcInletOutletBoundary, 0);
    addNamedToRunTimeSelectionTable
    (
        mcBoundary,
        mcInletOutletBoundary,
        mcBoundary,
        inletOutlet
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcInletOutletBoundary::mcInletOutletBoundary
(
    mcParticleCloud& cloud,
    label patchID,
    const dictionary& dict
)
:
    mcOpenBoundary(cloud, patchID, dict),
    probVec_(patch().size()),
    fwdTrans_(patch().size()),
    revTrans_(patch().size()),
    U_(),
    R_(),
    inRnd_()
{
    // compute triangle probability vector (normalized, cummulative sum of
    // areas), transformations and transform R_ into wall-normal system
    const polyPatch& pp = patch().patch();
    vector ex(1.0, 0.0, 0.0);
    const fvMesh& mesh = cloud.mesh();
    const vectorField& Sf = mesh.Sf().boundaryField()[patchID];
    const scalarField& magSf = mesh.magSf().boundaryField()[patchID];
    forAll(pp, faceI)
    {
        const face& f = pp[faceI];
        const point& C = mesh.C().boundaryField()[patchID][faceI];
        probVec_[faceI].setSize(f.size());

        pointField points = f.points(mesh.points());
        scalar area = 0;
        label n = points.size();
        for (label pointI1=0, pointI2=1; pointI1 < n; ++pointI1, ++pointI2)
        {
            if (pointI2==n) pointI2 = 0;
            point v1 = points[pointI1] - C;
            point v2 = points[pointI2] - C;
            area += mag(v1^v2); // don't bother with factor 0.5
            probVec_[faceI][pointI1] = area;
        }
        for (label pointI=0; pointI < n; ++pointI)
        {
            probVec_[faceI][pointI] /= area;
        }
        fwdTrans_[faceI] = rotationTensor(-Sf[faceI]/magSf[faceI], ex);
        revTrans_[faceI] = fwdTrans_[faceI].T();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField& Foam::mcInletOutletBoundary::getU()
{
    if (U_.empty())
    {
        const vectorField& U = cloud().Ufv().boundaryField()[patchID()];
        U_ = transform(fwdTrans_, U)();
    }
    return U_;
}


const Foam::symmTensorField& Foam::mcInletOutletBoundary::getR()
{
    if (R_.empty())
    {
        const symmTensorField& tau =
            cloud().TaucPdf().boundaryField()[patchID()];
        // clipping to compensate for bad divergence errors
        const compressible::turbulenceModel& tm = cloud().turbulenceModel();
        scalar k0;
        if (isA<compressible::RASModel>(tm))
        {
#if FOAM_HEX_VERSION < 0x200
            k0 = refCast<const compressible::RASModel>(tm).k0().value();
#else
            k0 = refCast<const compressible::RASModel>(tm).kMin().value();
#endif
        }
        else
        {
            k0 = SMALL;
        }
        R_ = max
        (
            transform(fwdTrans_, tau),
            symmTensor
            (
                k0, -GREAT, -GREAT,
                    k0,     -GREAT,
                            k0
            )
        );
    }
    return R_;
}


Foam::PtrList<Foam::mcInletRandom>&
Foam::mcInletOutletBoundary::getInRnd()
{
    if (inRnd_.empty())
    {
        const polyPatch& pp = patch().patch();
        inRnd_.setSize(pp.size());
        const vectorField& U = getU();
        const symmTensorField& R = getR();
        Random& rnd = cloud().random();
        forAll(pp, faceI)
        {
            inRnd_.set(faceI, mcInletRandom::New
                (
                    rnd,
                    U[faceI].x(),
                    sqrt(R[faceI].xx()),
                    subDict("randomGenerator")
                ));
        }
    }
    return inRnd_;
}


Foam::point Foam::mcInletOutletBoundary::randomPoint(label faceI)
{
    const fvMesh& mesh = cloud().mesh();
    const face& f = patch().patch()[faceI];
    const scalarList& w = probVec_[faceI];
    const point& C = mesh.C().boundaryField()[patchID()][faceI];
    const vector& Sf = mesh.Sf().boundaryField()[patchID()][faceI];
    pointField points = f.points(mesh.points());
    Random& rnd = cloud().random();
    // find triangle in which to place point
    label idx1 = std::distance(w.begin(), std::lower_bound(
            w.begin(), w.end(), rnd.scalar01()));
    label idx2 = idx1+1;
    if (idx2 == w.size())
    {
        idx2 = 0;
    }
    scalar a1 = rnd.scalar01();
    scalar a2 = rnd.scalar01();
    if (a1+a2 > 1.0)
    {
        // random point is outside triangle, mirror it inside
        a1 = 1.0 - a1;
        a2 = 1.0 - a2;
    }
    point p = a1*points[idx1] + a2*points[idx2] + (1.0-a1-a2)*C - SMALL*Sf;
    // If the case has reduced dimensionality, put the coordinate of the
    // reduced dimension onto the coordinate plane
    if (mesh.nGeometricD() <= 2)
    {
        meshTools::constrainDirection(mesh, mesh.geometricD(), p);
    }
    return p;
}


Foam::vector Foam::mcInletOutletBoundary::randomVelocity(label faceI)
{
    Random& rnd = cloud().random();
    const vector& U = getU()[faceI];
    const symmTensor& R = getR()[faceI];
    // Normal component has a special distributions, other components are Gaussian
    vector Up = vector
    (
        getInRnd()[faceI].value(),
        U.y() + sqrt(R.yy()) * rnd.GaussNormal(),
        U.z() + sqrt(R.zz()) * rnd.GaussNormal()
    );
    // Return Up in global coordinate system
    return transform(revTrans_[faceI], Up);
}


void Foam::mcInletOutletBoundary::correct(bool afterMove)
{
    mcOpenBoundary::correct(afterMove);

    if (afterMove)
    {
        return;
    }
    const fvMesh& mesh = cloud().mesh();
#ifdef FULLDEBUG
    if (debug > 1 && !statFile_.valid())
    {
        statFile_.reset
        (
            new OFstream
            (
                mesh.time().path() / patch().name()+".inletOutlet"
            )
        );
    }
#endif
    label np = 0;
    const polyPatch& pp = patch().patch();
    const scalarField& magSf = mesh.magSf().boundaryField()[patchID()];
    const scalarField& rho = cloud().rhocPdf().boundaryField()[patchID()];
    List<scalarField*> PhicPdf(cloud().PhicPdf().size());
    forAll(PhicPdf, PhiI)
    {
        PhicPdf[PhiI] = &cloud().PhicPdf()[PhiI]->boundaryField()[patchID()];
    }
    scalar dt = cloud().deltaT().value();
    const label Npc = cloud().solutionDict().particlesPerCell();
    const labelList& conservedScalars = cloud().conservedScalars();
    forAll(pp, faceI)
    {
        label cellI = pp.faceCells()[faceI];
        if (Un_[faceI] > 0)
        {
            // Mean velocity points out of the domain, so this is an outlet
            continue;
        }
        mcInletRandom& inrnd = getInRnd()[faceI];
        inrnd.updateCoeffs
        (
            getU()[faceI].x(),
            sqrt(getR()[faceI].xx())
        );
        // Mass of every particle generated
        scalar mp = rho[faceI]*mesh.V()[cellI]/Npc;
        // Mass flowing into domain across faceI during dt
        scalar mIn = rho[faceI]*inrnd.Q()*magSf[faceI]*dt;

        scalarField phi(PhicPdf.size());
        forAll(PhicPdf, PhiI)
        {
            phi[PhiI] = (*PhicPdf[PhiI])[faceI];
        }

        DynamicList<mcParticle*> genParticles;
        genParticles.reserve(mIn/mp);
        scalar mGen = 0;
        while (mGen < mIn)
        {
            vector Up = randomVelocity(faceI);
            mcParticle* p = new mcParticle
            (
                cloud(),
                randomPoint(faceI),
                cellI,
                mp,
                Up,
                phi
            );
            p->isOnInletBoundary() = true;
            cloud().localTimeStepping().correct(*p);
            p->m() /= p->eta();
            if (mGen + p->m() > mIn)
            {
                if (cloud().random().scalar01() > (mIn - mGen)/p->m())
                {
                    delete p;
                    break;
                }
            }
            cloud().addParticle(p);
            genParticles.append(p);
            mGen += p->m();
            ++np;
#ifdef FULLDEBUG
            if (debug > 1)
            {
                statFile_() << Up.x() << tab << Up.y() << tab << Up.z()
                    << tab << mp << nl;
            }
#endif
        }
        if (genParticles.size())
        {
            cloud().adjustAxiSymmetricMass(genParticles);
            forAll(genParticles, i)
            {
                mcParticle& p = *genParticles[i];
                scalar mpd = cloud().massPerDepth(p);
                massIn()[0] += mpd;
                forAll(conservedScalars, csI)
                {
                    massIn()[csI+1] += mpd*p.Phi()[conservedScalars[csI]];
                }
            }
        }
    }
    if (debug)
    {
        Pout<< "Inlet particles: " << np << " generated"
            " on patch " << patch().name() << endl;
    }
#ifdef FULLDEBUG
    if (debug > 1)
    {
        statFile_().flush();
    }
#endif
}

// ************************************************************************* //
