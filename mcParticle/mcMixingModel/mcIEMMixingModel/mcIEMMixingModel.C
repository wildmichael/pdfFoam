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

#include "mcIEMMixingModel.H"

#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcIEMMixingModel, 0);
    addNamedToRunTimeSelectionTable
    (
        mcMixingModel,
        mcIEMMixingModel,
        mcMixingModel,
        IEM
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcIEMMixingModel::mcIEMMixingModel
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& dict
)
:
    mcMixingModel(db, dict),
    Cmix2_(0.5*lookupOrDefault<scalar>("Cmix", 2.0))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcIEMMixingModel::correct(Foam::mcParticleCloud& cloud)
{
    const labelList& mixedScalars = cloud.mixedScalars();
    label nMix = mixedScalars.size();
    // prepare interpolators for mean fields
    PtrList<interpolationCellPoint<scalar> > PhiInterp(nMix);
    forAll(mixedScalars, mixedI)
    {
        interpolationCellPoint<scalar>* interp =
            new interpolationCellPoint<scalar>
            (
                *cloud.PhicPdf()[mixedScalars[mixedI]]
            );
        PhiInterp.set(mixedI, interp);
    }
    const fvMesh& mesh = cloud.mesh();
    scalar dt =  mesh.time().deltaT().value();
    // loop over particles
    forAllIter(mcParticleCloud, cloud, pIter)
    {
        mcParticle& p = pIter();
        cellPointWeight cpw(mesh, p.position(), p.cell(), p.face());
        scalar Omega = p.Omega();
        // loop over mixed properties
        forAll(mixedScalars, mixedI)
        {
            // interpolate property to particle location
            scalar phiap = PhiInterp[mixedI].interpolate(cpw);
            // apply IEM mixing
            scalar& phi = p.Phi()[mixedScalars[mixedI]];
            phi -= Cmix2_ * Omega * (phi - phiap) * dt;
        }
    }
}

// ************************************************************************* //
