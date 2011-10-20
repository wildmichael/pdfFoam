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
    scalar Umean,
    scalar urms,
    const dictionary& dict
)
:
    mcInletRandom(rnd, Umean, urms, dict)
{
    updateCoeffs(Umean, urms);
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
    b_ = M_SQRT1_2/urms();
    b2_ = b_*b_;
    expma2b2_ = exp(-Umean()*Umean()*b2_);
    abSqrtPi_ = Umean()*b_*sqrt(M_PI);
    erfab_ = erf(Umean()*b_);
    denom_ = expma2b2_+abSqrtPi_*(1.0+erfab_);
    b22denom_ = 2.0*b2_/denom_;
    m_ = 0.5*(Umean()*b_+sqrt(2.0+Umean()*Umean()*b2_))/b_;
}

// ************************************************************************* //
