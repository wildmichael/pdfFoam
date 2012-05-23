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

#include "mcNaiveFlameletModel.H"

#include "addToRunTimeSelectionTable.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcNaiveFlameletModel, 0);
    addNamedToRunTimeSelectionTable
    (
        mcReactionModel,
        mcNaiveFlameletModel,
        mcReactionModel,
        naiveFlamelet
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcNaiveFlameletModel::mcNaiveFlameletModel
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcReactionModel(cloud, db, subDictName),
    zIdx_(findIdx("zName", "z")),
    TIdx_(findIdx("TName", "T"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcNaiveFlameletModel::correct(mcParticle& p)
{
    const scalar& z = p.Phi()[zIdx_];
    p.rho() = (1. - 3.2*z*(1.-z));
    p.Phi()[TIdx_] = 1e5/p.rho();
}

// ************************************************************************* //
