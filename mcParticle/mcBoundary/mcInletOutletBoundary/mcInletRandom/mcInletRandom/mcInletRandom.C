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
    scalar UMean,
    scalar uRms,
    const dictionary& dict
)
:
    dictionary(dict),
    rnd_(rnd)
{
    updateCoeffs(UMean, uRms);
}

// * * * * * * * * * * * * * * * * Destructor *  * * * * * * * * * * * * * * //

Foam::mcInletRandom::~mcInletRandom()
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mcInletRandom> Foam::mcInletRandom::New
(
    Random& rnd,
    scalar UMean,
    scalar uRms,
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
    return autoPtr<mcInletRandom>(cstrIter()(rnd, UMean, uRms, dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcInletRandom::updateCoeffs
(
    Foam::scalar UMean,
    Foam::scalar uRms
)
{
    UMean_ = UMean;
    uRms_ = uRms;
    b_ = M_SQRT1_2/uRms_;
    b2_ = b_*b_;
    expma2b2_ = exp(-UMean_*UMean_*b2_);
    abSqrtPi_ = UMean_*b_*sqrt(M_PI);
    erfab_ = erf(UMean_*b_);
    denom_ = expma2b2_+abSqrtPi_*(1.+erfab_);
    b22denom_ = 2.*b2_/denom_;
    UCondMean_ = 0.5*(UMean_+expma2b2_/(b_*sqrt(M_PI))+UMean_*erfab_);
}

// ************************************************************************* //
