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

#include "mcInversionInletRandom.H"
#include "addToRunTimeSelectionTable.H"
#include "Random.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcInversionInletRandom, 0);
    addNamedToRunTimeSelectionTable
    (
        mcInletRandom,
        mcInversionInletRandom,
        mcInletRandom,
        inversion
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcInversionInletRandom::mcInversionInletRandom
(
    Random& rnd,
    scalar UMean,
    scalar uRms,
    const dictionary& dict
)
:
    mcInletRandom(rnd, UMean, uRms, dict)
{
    updateCoeffs(UMean, uRms);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::mcInversionInletRandom::newton
(
    scalar x0,
    scalar y,
    label& nIter
)
{
    static const label nMaxIter = 50;
    static const scalar tol = 1e-8;
    nIter = 0;
    do
    {
      scalar fval = CDF(x0) - y;
      scalar fder = PDF(x0);
      if (fabs(fder) < VSMALL)
      {
          WarningIn("mcInversionInletRandom::newton(scalar, scalar, label&)")
              << "Derivative was zero (" << fder << ") after "
              << nIter << " iterations."
              << endl;
        return x0;
      }
      scalar x = x0 - fval/fder;
      if (fabs(x - x0) < tol)
      {
          return x;
      }
      x0 = x;
      ++nIter;
    } while (nIter < nMaxIter);
    WarningIn("mcInversionInletRandom::newton(scalar, scalar, label&)")
        << "Exceeded maxinum number of iterations (" << nMaxIter << ")."
        << endl;
    return x0;
}


Foam::scalar Foam::mcInversionInletRandom::value()
{
    scalar U = rnd().scalar01();
    // perturbation away from peak position
    scalar d = sign(U-CDF(m_))*SMALL;
    label nIter;
    scalar r = newton(m_+d, U, nIter);
    return r;
}


void Foam::mcInversionInletRandom::updateCoeffs
(
    scalar U,
    scalar up
)
{
    mcInletRandom::updateCoeffs(U, up);
    m_ = 0.5*(UMean()*b()+sqrt(2.0+UMean()*UMean()*b2()))/b();
}

// ************************************************************************* //
