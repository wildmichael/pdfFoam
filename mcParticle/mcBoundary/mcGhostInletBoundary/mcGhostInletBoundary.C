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

#include "mcGhostInletBoundary.H"

#include "addToRunTimeSelectionTable.H"
#include "mcParticle.H"
#include "mcParticleCloud.H"
#include "surfaceMesh.H"
#include "fvsPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcGhostInletBoundary, 0);
    addNamedToRunTimeSelectionTable
    (
        mcBoundary,
        mcGhostInletBoundary,
        mcBoundary,
        ghostInlet
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcGhostInletBoundary::mcGhostInletBoundary
(
    const Foam::fvMesh& mesh,
    Foam::label patchID,
    const Foam::dictionary& dict
)
:
    mcGhostBoundary(mesh, patchID, dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcGhostInletBoundary::correct
(
    Foam::mcParticleCloud& cloud,
    bool afterMove
)
{
    mcGhostBoundary::correct(cloud, afterMove);
}

void Foam::mcGhostInletBoundary::hitPatch
(
    Foam::mcParticle& p,
    Foam::mcParticle::trackData& td
)
{
    // TODO moving-wall reflection
    // (Meyer and Jenny 2007: http://dx.doi.org/10.1016/j.jcp.2007.06.014)
    td.keepParticle = false;
}


void Foam::mcGhostInletBoundary::hitPatch
(
    Foam::mcParticle& p,
    int&
)
{}

// ************************************************************************* //
