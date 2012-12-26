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

#include "mcSamplingInletRandom.H"

#include "addToRunTimeSelectionTable.H"
#include "Random.H"
#include "scalarList.H"

namespace {

// * * * * * * * * * * * * * Local Helper Functions * * * * * * * * * * * * //

using namespace Foam;
void fillBuf
(
    Random&            rnd,
    scalar             UMean,
    scalar             uRms,
    FIFOStack<scalar>& buf
)
{
    static scalarList tmp(1e6);
    static const scalar CFL = 0.1;
    scalar m = -VGREAT;
    forAll(tmp, i)
    {
        tmp[i] = rnd.GaussNormal()*uRms + UMean;
        m = max(m, fabs(tmp[i]));
    }
    forAll(tmp, i)
    {
        if (rnd.scalar01() + tmp[i]*CFL/m > 1.)
        {
            buf.push(tmp[i]);
        }
    }
}

} // End namespace (anonymous)

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcSamplingInletRandom, 0);
    addNamedToRunTimeSelectionTable
    (
        mcInletRandom,
        mcSamplingInletRandom,
        mcInletRandom,
        sampling
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcSamplingInletRandom::mcSamplingInletRandom
(
    Random& rnd,
    scalar UMean,
    scalar uRms,
    const dictionary& dict
)
:
    mcInletRandom(rnd, UMean, uRms, dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::mcSamplingInletRandom::value()
{
    if (!buf_.size())
    {
        fillBuf(rnd(), UMean(), uRms(), buf_);
    }
    return buf_.pop();
}


void Foam::mcSamplingInletRandom::updateCoeffs
(
    scalar U,
    scalar up
)
{
    if (fabs(UMean()-U)>VSMALL || fabs(uRms()-up)>VSMALL)
    {
        buf_.clear();
    }
    mcInletRandom::updateCoeffs(U, up);
}

// ************************************************************************* //
