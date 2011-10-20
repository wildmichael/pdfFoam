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

#include "mcInletRandom.H"

#include "Random.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcInletRandom, 0);
    defineRunTimeSelectionTable(mcInletRandom, mcInletRandom);

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcInletRandom::mcInletRandom
(
    Foam::Random& rnd,
    scalar Umean,
    scalar urms,
    const dictionary& dict
)
:
    dictionary(dict),
    rnd_(rnd)
{
    updateCoeffs(Umean, urms);
}

// * * * * * * * * * * * * * * * * Destructor *  * * * * * * * * * * * * * * //

Foam::mcInletRandom::~mcInletRandom()
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mcInletRandom> Foam::mcInletRandom::New
(
    Random& rnd,
    scalar Umean,
    scalar urms,
    const dictionary& dict
)
{
    word randomType(dict.lookup("type"));

    mcInletRandomConstructorTable::iterator cstrIter =
        mcInletRandomConstructorTablePtr_->find(randomType);

    if (cstrIter == mcInletRandomConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mcInletRandom::New(Random&, scalar, scalar, const dictionary&)"
        )   << "Unknown mcInletRandom type " << randomType << endl << endl
            << "Valid mcInletRandom types are :" << endl
            << mcInletRandomConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<mcInletRandom>(cstrIter()(rnd, Umean, urms, dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcInletRandom::updateCoeffs
(
    Foam::scalar Umean,
    Foam::scalar urms
)
{
    Umean_ = Umean;
    urms_ = urms;
}

// ************************************************************************* //
