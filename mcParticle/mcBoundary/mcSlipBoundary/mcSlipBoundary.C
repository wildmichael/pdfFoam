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

#include "mcSlipBoundary.H"

#include "addToRunTimeSelectionTable.H"
#include "mcParticle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcSlipBoundary, 0);
    addNamedToRunTimeSelectionTable
    (
        mcBoundary,
        mcSlipBoundary,
        mcBoundary,
        slip
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcSlipBoundary::mcSlipBoundary
(
    const Foam::fvMesh& mesh,
    Foam::label patchID,
    const Foam::dictionary& dict
)
:
    mcBoundary(mesh, patchID, dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcSlipBoundary::hitPatch
(
    Foam::mcParticle& p,
    Foam::mcParticle::trackData& td
)
{
    vector nw = patch().faceAreas()[patch().whichFace(p.face())];
    nw /= mag(nw); // Wall normal (outward)
    p.transformProperties(I - 2.0*nw*nw);
    p.reflected() = true;
}


void Foam::mcSlipBoundary::hitPatch
(
    Foam::mcParticle& p,
    int&
)
{}

// ************************************************************************* //
