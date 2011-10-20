/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mcInletOutletBoundary.H"

#include "addToRunTimeSelectionTable.H"
#include "mcParticle.H"
#include "mcParticleCloud.H"
#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "Random.H"
#include "transformField.H"
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
    const Foam::fvMesh& mesh,
    Foam::label patchID,
    const Foam::dictionary& dict
)
:
    mcOpenBoundary(mesh, patchID, dict),
    probVec_(patch().size()),
    fwdTrans_(patch().size()),
    revTrans_(patch().size()),
    U_(),
    urms_(),
    inRnd_()
{
    // compute triangle probability vector (normalized, cummulative sum of
    // areas), transformations and transform urms_ into wall-normal system
    const polyPatch& pp = patch();
    vector ex(1.0, 0.0, 0.0);
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
        fwdTrans_[faceI] = rotationTensor(ex, -Sf[faceI]/magSf[faceI]);
        revTrans_[faceI] = fwdTrans_[faceI].T();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField& Foam::mcInletOutletBoundary::getU
(
    const Foam::mcParticleCloud& cloud
)
{
    if (U_.empty())
    {
        const vectorField& U = cloud.Ufv().boundaryField()[patchID()];
        U_ = transform(fwdTrans_, U)();
    }
    return U_;
}


const Foam::vectorField& Foam::mcInletOutletBoundary::geturms
(
    const Foam::mcParticleCloud& cloud
)
{
    if (urms_.empty())
    {
        const symmTensorField& tau = cloud.TaucPdf().boundaryField()[patchID()];
        vectorField urms = vector(1., 0., 0.)*sqrt(tau.component(symmTensor::XX))
                         + vector(0., 1., 0.)*sqrt(tau.component(symmTensor::YY))
                         + vector(0., 0., 1.)*sqrt(tau.component(symmTensor::ZZ));
        urms_ = transform(fwdTrans_, urms)();
    }
    return urms_;
}


Foam::PtrList<Foam::InletRandom>&
Foam::mcInletOutletBoundary::getInRnd
(
    Foam::mcParticleCloud& cloud
)
{
    if (inRnd_.empty())
    {
        const polyPatch& pp = patch();
        inRnd_.setSize(pp.size());
        const vectorField& U = getU(cloud);
        const vectorField& urms = geturms(cloud);
        Random& rnd = cloud.random();
        forAll(pp, faceI)
        {
            inRnd_.set(faceI, new InletRandom
                (
                    rnd,
                    U[faceI].x(),
                    urms[faceI].x()
                ));
        }
    }
    return inRnd_;
}


Foam::point Foam::mcInletOutletBoundary::randomPoint
(
    Foam::mcParticleCloud& cloud,
    Foam::label faceI
)
{
    const face& f = patch()[faceI];
    const scalarList& w = probVec_[faceI];
    const point& C = mesh().C().boundaryField()[patchID()][faceI];
    const vector& Sf = mesh().Sf().boundaryField()[patchID()][faceI];
    pointField points = f.points(mesh().points());
    Random& rnd = cloud.random();
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
    return a1*points[idx1] + a2*points[idx2] + (1.0-a1-a2)*C - SMALL*Sf;
}


Foam::vector Foam::mcInletOutletBoundary::randomVelocity
(
    Foam::mcParticleCloud& cloud,
    Foam::label faceI
)
{
    Random& rnd = cloud.random();
    const vector& U = getU(cloud)[faceI];
    const vector& urms = geturms(cloud)[faceI];
    // Normal component has a special distributions, other components are Gaussian
    vector Up = vector
    (
        getInRnd(cloud)[faceI].value(),
        U.y() + urms.y() * rnd.GaussNormal(),
        U.z() + urms.z() * rnd.GaussNormal()
    );
    // Return Up in global coordinate system
    return transform(revTrans_[faceI], Up);
}


void Foam::mcInletOutletBoundary::correct
(
    Foam::mcParticleCloud& cloud,
    bool afterMove
)
{
    if (afterMove)
    {
        return;
    }
#ifdef FULLDEBUG
    if (debug > 1 && !statFile_.valid())
    {
        statFile_.reset
        (
            new OFstream
            (
                mesh().time().path() / patch().name()+".inletOutlet"
            )
        );
    }
#endif
    label np = 0;
    const polyPatch& pp = patch();
    const vectorField& Sf = mesh().Sf().boundaryField()[patchID()];
    const vectorField& U = cloud.Ufv().boundaryField()[patchID()];
    const scalarField& rho = cloud.rhocPdf().boundaryField()[patchID()];
    List<scalarField*> PhicPdf(cloud.PhicPdf().size());
    forAll(PhicPdf, PhiI)
    {
        PhicPdf[PhiI] = &cloud.PhicPdf()[PhiI]->boundaryField()[patchID()];
    }
    Random& rnd = cloud.random();
    scalar dt = mesh().time().deltaT().value();
    scalar Npc = cloud.Npc();
    forAll(pp, faceI)
    {
        label cellI = pp.faceCells()[faceI];
        scalar vol = -Sf[faceI] & U[faceI] * dt;
        if (vol < 0)
        {
            // Mean velocity points out of the domain, so this is an outlet
            continue;
        }
        getInRnd(cloud)[faceI].updateCoeffs(getU(cloud)[faceI].x(),
                                            geturms(cloud)[faceI].x());
        // number of particles to generate
        scalar Nf = vol / mesh().V()[cellI] * Npc;
        label N = floor(Nf);
        N += rnd.scalar01() < (Nf-N);
        if (N < 1)
        {
            continue;
        }
        np += N;
        // mass of every particle
        scalar m = vol * rho[faceI] / N;

        scalarField phi(PhicPdf.size());
        forAll(PhicPdf, PhiI)
        {
            phi[PhiI] = (*PhicPdf[PhiI])[faceI];
        }

        for (label i=0; i<N; ++i)
        {
            vector u = randomVelocity(cloud, faceI);
            mcParticle* p = new mcParticle
            (
                cloud,
                randomPoint(cloud, faceI),
                cellI,
                m,
                U[faceI],
                u,
                U[faceI],
                phi
            );
            p->isOnInletBoundary() = true;
            cloud.addParticle(p);
#ifdef FULLDEBUG
            if (debug > 1)
            {
                statFile_() << u.x() << tab << u.y() << tab << u.z() << nl;
            }
#endif
        }
#ifdef FULLDEBUG
            if (debug > 1)
            {
                statFile_().flush();
            }
#endif
        if (debug)
        {
            reduce(np, sumOp<label>());
            Info<< "Inlet particles: " << np << " generated"
                " on patch " << patch().name() << endl;
        }
    }
}

// ************************************************************************* //
