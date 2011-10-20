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

#include "mcBoundary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcBoundary, 0);
    defineRunTimeSelectionTable(mcBoundary, mcBoundary);

} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcBoundary::mcBoundary
(
    const Foam::fvMesh& mesh,
    label patchID,
    const Foam::dictionary& dict
)
:
    dictionary(dict),
    mesh_(mesh),
    patchID_(patchID),
    patch_(mesh.boundaryMesh()[patchID])
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mcBoundary> Foam::mcBoundary::New
(
    const Foam::fvMesh& mesh,
    label patchID,
    const Foam::dictionary& dict
)
{
    word boundaryType(dict.lookup("type"));

    mcBoundaryConstructorTable::iterator cstrIter =
        mcBoundaryConstructorTablePtr_->find(boundaryType);

    if (cstrIter == mcBoundaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mcBoundary::New(const fvMesh&, label, const dictionary&)"
        )   << "Unknown mcBoundary type " << boundaryType << endl << endl
            << "Valid mcBoundary types are :" << endl
            << mcBoundaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<mcBoundary>(cstrIter()(mesh, patchID, dict));
}

// ************************************************************************* //
