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

#include "mcParticleLocalTimeStepping.H"

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "interpolationCellPointFace.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcParticleLocalTimeStepping, 0);
    addNamedToRunTimeSelectionTable
    (
        mcLocalTimeStepping,
        mcParticleLocalTimeStepping,
        mcLocalTimeStepping,
        particle
    );

} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcParticleLocalTimeStepping::mcParticleLocalTimeStepping
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& dict
)
:
    mcLocalTimeStepping(db, dict),
    CourantU_
    (
        max(lookupOrAddDefault<scalar>("CourantU", 0.3), 1e-6)
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcParticleLocalTimeStepping::correct
(
    Foam::mcParticleCloud& cloud
)
{
    forAllIter(mcParticleCloud, cloud, pIter)
    {
        mcParticle& p = pIter();
        p.eta() = p.Co() > SMALL ? CourantU_ / p.Co() : 1.0;
    }
}

// ************************************************************************* //
