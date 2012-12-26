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

#include "mcOpenBoundary.H"

#include "mcParticle.H"
#include "mcParticleCloud.H"
#include "surfaceMesh.H"
#include "fvsPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcOpenBoundary, 0);

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcOpenBoundary::mcOpenBoundary
(
    mcParticleCloud& cloud,
    label patchID,
    const dictionary& dict
)
:
    mcBoundary(cloud, patchID, dict),
    reflecting_(lookupOrDefault<Switch>("reflecting", true)),
    n_
    (
        cloud.mesh().Sf().boundaryField()[patchID]
      / cloud.mesh().magSf().boundaryField()[patchID]
    ),
    Un_(patch().size(), 0.)
{
    const polyPatch& pp = patch().patch();
    forAll(pp, faceI)
    {
        label cellI = pp.faceCells()[faceI];
        if (!cellFaces_.found(cellI))
        {
            cellFaces_.insert(cellI, new DynamicList<label>());
        }
        cellFaces_[cellI]->append(faceI);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcOpenBoundary::correct(bool afterMove)
{
    if (!reflecting_)
    {
        return;
    }

    // Guard to ensure that the reflected particles are only treated by one
    // instance of mcOpenBoundary globally on one processor.
    static bool handledReflectedParticles = false;
    if (!afterMove)
    {
        const fvMesh& mesh = cloud().mesh();
        handledReflectedParticles = false;
        const label patchI = patchID();
        const vectorField& Cf     = mesh.Cf().boundaryField()[patchI];
        const scalar dt           = cloud().deltaT().value();
        const labelList& conservedScalars = cloud().conservedScalars();
        Un_ =
            (
                n_
              & patch().lookupPatchField<volVectorField, vector>
                (
                    lookupOrDefault<word>("UName", "U")
                )
            );

        // Estimate location of "moving boundary" at t=t0
        scalarList x0 = Un_*dt;

        // Discard particles behind the "moving boundary"
        forAllIter(mcParticleCloud, cloud(), pIter)
        {
            mcParticle& p = pIter();
            label cellI = p.cell();
            dynamicLabelListPtrMap::const_iterator cfIter
            (
                cellFaces_.find(cellI)
            );
            // Is particle in cell adjacent to this patch?
            if (cfIter != cellFaces_.end())
            {
                const DynamicList<label>& faces = *cfIter();
                forAll(faces, i)
                {
                    label faceI = faces[i];
                    // outflow? (intentionally "not inflow")
                    if (!(Un_[faceI] < 0.))
                    {
                        scalar xp = n_[faceI]&(Cf[faceI]-p.position());
                        if (xp < p.eta()*x0[faceI])
                        {
                            scalar mpd = cloud().massPerDepth(p);
                            massOut()[0] -= mpd;
                            forAll(conservedScalars, csI)
                            {
                                label i = conservedScalars[csI];
                                massOut()[csI+1] -= mpd*p.Phi()[i];
                            }
                            // FIXME Does deletion invalidate the iteration?
                            // Not sure how IDLListBase implements this...
                            cloud().deleteParticle(p);
                        }
                    }
                }
            }
        }
    }
    else if (!handledReflectedParticles)
    {
        // Update velocities of all reflected particles
        forAllIter(mcParticleCloud, cloud(), pIter)
        {
            mcParticle& p = pIter();
            if (p.reflectedAtOpenBoundary())
            {
                p.UParticle() += 2*p.reflectionBoundaryVelocity();
                p.reflectedAtOpenBoundary() = false;
            }
        }
        handledReflectedParticles = true;
    }
}


void Foam::mcOpenBoundary::hitPatch
(
    mcParticle& p,
#if FOAM_HEX_VERSION < 0x200
    mcParticle::trackData& td
#else
    mcParticle::trackData& td,
    const label patchI,
    const scalar trackFraction,
    const tetIndices& tetIs
#endif
)
{
    const polyPatch& pp = patch().patch();
    label faceI = pp.whichFace(p.face());
    //const labelList& conservedScalars = cloud().conservedScalars();
    //
    // If this boundary is not reflecting or if it is an inflow boundary,
    // delete particle.
    //if (Un_[faceI] < 0)
    //{
    //    scalar mpd = cloud().massPerDepth(p);
    //    massIn()[0] -= mpd;
    //    forAll(conservedScalars, csI)
    //    {
    //        massIn()[csI+1] -= mpd*p.Phi()[conservedScalars[csI]];
    //    }
    //    td.keepParticle = false;
    //    return;
    //}
    //if (!reflecting_)
    //{
    //    scalar mpd = cloud().massPerDepth(p);
    //    massOut()[0] -= mpd;
    //    forAll(conservedScalars, csI)
    //    {
    //        massOut()[csI+1] -= mpd*p.Phi()[conservedScalars[csI]];
    //    }
    //    td.keepParticle = false;
    //    return;
    //}

    p.transformProperties(I - 2.0*n_[faceI]*n_[faceI]);
    p.reflected() = true;
    p.reflectedAtOpenBoundary() = true;
    p.reflectionBoundaryVelocity() = n_[faceI]*Un_[faceI];
}


#if FOAM_HEX_VERSION < 0x200
void Foam::mcOpenBoundary::hitPatch
(
    Foam::mcParticle& p,
    int&
)
{}
#endif

// ************************************************************************* //
