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
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcLocalTimeStepping(cloud, db, subDictName),
    CourantU_(0.),
    upperBound_(0.)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcParticleLocalTimeStepping::updateInternals()
{
    CourantU_ =
        max(solutionDict().lookupOrDefault<scalar>("CourantU", 0.3), 1e-6);
    upperBound_ =
        solutionDict().lookupOrDefault<scalar>("upperBound", 10);
}


void Foam::mcParticleLocalTimeStepping::correct(mcParticle& p)
{
    p.eta() = min(CourantU_ / stabilise(p.Co(), VSMALL), upperBound_);
}

// ************************************************************************* //
