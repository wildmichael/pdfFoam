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

#include "mcMixingModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcMixingModel, 0);
    defineRunTimeSelectionTable(mcMixingModel, mcMixingModel);

} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcMixingModel::mcMixingModel
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& dict
)
:
    mcModel(db, dict)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mcMixingModel> Foam::mcMixingModel::New
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& dict
)
{
    word mixingType(dict.lookup("mixingModel"));

    mcMixingModelConstructorTable::iterator cstrIter =
        mcMixingModelConstructorTablePtr_->find(mixingType);

    if (cstrIter == mcMixingModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mcMixingModel::New(const fvMesh&, const dictionary&)"
        )   << "Unknown mcMixingModel type " << mixingType << endl << endl
            << "Valid mcMixingModel types are :" << endl
            << mcMixingModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    const dictionary& coeffs = dict.subDict(mixingType+"Coeffs");
    return autoPtr<mcMixingModel>(cstrIter()(db, coeffs));
}

// ************************************************************************* //
