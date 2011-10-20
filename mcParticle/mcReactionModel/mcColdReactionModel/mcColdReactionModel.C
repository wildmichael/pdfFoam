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

#include "mcColdReactionModel.H"

#include "addToRunTimeSelectionTable.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcColdReactionModel, 0);
    addNamedToRunTimeSelectionTable
    (
        mcReactionModel,
        mcColdReactionModel,
        mcReactionModel,
        cold
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcColdReactionModel::mcColdReactionModel
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& dict
)
:
    mcReactionModel(db, dict),
    rho_(dict.lookupOrDefault<scalar>("density", 1.0))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcColdReactionModel::correct(Foam::mcParticleCloud& cloud)
{
    forAllIter(mcParticleCloud, cloud, pIter)
    {
        correct(cloud, pIter());
    }
}


void Foam::mcColdReactionModel::correct
(
    Foam::mcParticleCloud& cloud,
    Foam::mcParticle& p
)
{
    p.rho() = rho_;
}

// ************************************************************************* //
