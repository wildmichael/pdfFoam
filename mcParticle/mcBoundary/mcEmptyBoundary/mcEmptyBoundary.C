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

#include "mcEmptyBoundary.H"

#include "addToRunTimeSelectionTable.H"
#include "mcParticle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcEmptyBoundary, 0);
    addNamedToRunTimeSelectionTable
    (
        mcBoundary,
        mcEmptyBoundary,
        mcBoundary,
        empty
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcEmptyBoundary::mcEmptyBoundary
(
    const Foam::fvMesh& mesh,
    Foam::label patchID,
    const Foam::dictionary& dict
)
:
    mcBoundary(mesh, patchID, dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcEmptyBoundary::hitPatch
(
    Foam::mcParticle& p,
    Foam::mcParticle::trackData& td
)
{
    FatalErrorIn("mcEmptyBoundary::hitPatch"
        "(mcParticle&, mcParticle::trackData&)")
        << "A mcParticle must never hit an empty boundary!" << endl
        << Foam::exit(FatalError);
}


void Foam::mcEmptyBoundary::hitPatch
(
    Foam::mcParticle& p,
    int&
)
{
    FatalErrorIn("mcEmptyBoundary::hitPatch"
        "(mcParticle&, int&)")
        << "A mcParticle must never hit an empty boundary!" << endl
        << Foam::exit(FatalError);
}

// ************************************************************************* //
