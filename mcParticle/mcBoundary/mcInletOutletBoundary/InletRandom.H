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

Class
    Foam::InletRandom

Description
    Generate wall-normal particle-velocity for an inlet boundary

    Generate random values from the PDF
    \f[ f_x = \frac{X}{N} e^{-(X-a)^2 b^2} \f]
    where
    \f{eqnarray*}{
        b &=& \frac{1}{\sqrt{2}u_\mathrm{rms}} \quad \text{and} \\
        N &=& b^2 e^{-\frac{a^2}{2 b^2}} + a \sqrt{\frac{\pi b^2}{2}}
              \left(1+\mathrm{erf}\left(\frac{a}{\sqrt{2 b^2}}\right)\right)
    \f}
    using numerical inversion of the CDF with Newton-Raphson.


SourceFiles
    InletRandom.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef InletRandom_H
#define InletRandom_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Random;

/*---------------------------------------------------------------------------*\
                         Class InletRandom Declaration
\*---------------------------------------------------------------------------*/

class InletRandom
{
    // Private data

        //- Random number generator
        Random& rnd_;
        //- Mean
        scalar a_;
        //- 1/sqrt(2*Variance)
        scalar b_;
        //- Commonly used, precomputed factors
        scalar b2_, expma2b2_, abSqrtPi_, erfab_, denom_, b22denom_;
        //- Location of the peak of the PDF
        scalar m_;

    // Private member functions

        //- Newton-Raphson root finding of the CDF
        scalar newton(scalar x0, scalar y, label& nIter);

        // Disallow default bitwise copy construct and assignment
        InletRandom(const InletRandom&);
        void operator=(const InletRandom&);

public:

    // Constructors

        //- Construct a random number generator given a Random object and
        //the mean velocity and the RMS of the fluctuations
        InletRandom(Random& rnd, scalar Umean, scalar urms);

    // Public member functions

        //- Generate one random number
        scalar value();

        //- Update the coefficients for changed Umean and urms
        void updateCoeffs(scalar Umean, scalar urms);

        //- Cummulative distribution function (CDF)
        inline scalar CDF(scalar x);

        //- Probability density function (PDF)
        inline scalar PDF(scalar x);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InletRandomI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
