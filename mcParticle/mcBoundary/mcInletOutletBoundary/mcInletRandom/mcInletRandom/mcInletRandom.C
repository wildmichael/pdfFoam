/*---------------------------------------------------------------------------*\
                pdfFoam: General Purpose PDF Solution Algorithm
                   for Reactive Flow Simulations in OpenFOAM

 Copyright (C) 2012 Michael Wild, Heng Xiao, Patrick Jenny,
                    Institute of Fluid Dynamics, ETH Zurich
-------------------------------------------------------------------------------
License
    This file is part of pdfFoam.

    pdfFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) version 3 of the same License.

    pdfFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with pdfFoam.  If not, see <http://www.gnu.org/licenses/>.

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
    expmU2b2_ = exp(-UMean_*UMean_*b2_);
    UbSqrtPi_ = UMean_*b_*sqrt(M_PI);
    erfUb_ = erf(UMean_*b_);
    denom_ = expmU2b2_+UbSqrtPi_*(1.+erfUb_);
    b22denom_ = 2.*b2_/denom_;
    Q_ = 0.5*(UMean_ + expmU2b2_/(b_*sqrt(M_PI)) + UMean_*erfUb_);
}

// ************************************************************************* //
