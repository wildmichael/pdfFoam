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
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcMixingModel(cloud, db, subDictName),
    Cmix2_(0.)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcIEMMixingModel::updateInternals()
{
    Cmix2_ = 0.5*solutionDict().lookupOrDefault<scalar>("Cmix", 2.0);
}


void Foam::mcIEMMixingModel::correct(mcParticle& p)
{
    const labelList& mixedScalars = cloud().mixedScalars();
    const scalar& dt =  cloud().deltaT().value();
    const scalar& Omega = p.Omega();
    const scalar& eta = p.eta();
    label cellI = p.cell();
    // loop over mixed properties
    forAll(mixedScalars, mixedI)
    {
        // apply IEM mixing
        scalar& phi = p.Phi()[mixedScalars[mixedI]];
        scalar& phiap = (*cloud().PhicPdf()[mixedScalars[mixedI]])[cellI];
        phi -= (1.0 - exp(-Cmix2_*Omega*eta*dt))*(phi - phiap);
    }
}

// ************************************************************************* //
