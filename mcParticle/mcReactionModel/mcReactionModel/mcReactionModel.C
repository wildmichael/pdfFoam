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

#include "mcReactionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcReactionModel, 0);
    defineRunTimeSelectionTable(mcReactionModel, mcReactionModel);

} // namespace Foamm

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcReactionModel::mcReactionModel
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& dict
)
:
    mcModel(db, dict)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mcReactionModel> Foam::mcReactionModel::New
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& dict
)
{
    word reactionType(dict.lookup("reactionModel"));

    mcReactionModelConstructorTable::iterator cstrIter =
        mcReactionModelConstructorTablePtr_->find(reactionType);

    if (cstrIter == mcReactionModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mcReactionModel::New(const fvMesh&, const dictionary&)"
        )   << "Unknown mcReactionModel type " << reactionType << endl << endl
            << "Valid mcReactionModel types are :" << endl
            << mcReactionModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<mcReactionModel>(cstrIter()(db, dict));
}

// ************************************************************************* //
