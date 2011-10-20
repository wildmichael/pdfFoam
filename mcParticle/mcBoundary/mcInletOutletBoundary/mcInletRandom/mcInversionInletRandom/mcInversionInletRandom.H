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
    Foam::mcInversionInletRandom

Description
    Generate wall-normal particle-velocity for an inlet boundary using
    numerical inversion of the CDF with Newton-Raphson.


SourceFiles
    mcInversionInletRandom.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcInversionInletRandom_H
#define mcInversionInletRandom_H

#include "mcInletRandom.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Random;

/*---------------------------------------------------------------------------*\
                    Class mcInversionInletRandom Declaration
\*---------------------------------------------------------------------------*/

class mcInversionInletRandom
    : public mcInletRandom
{
    // Private data

        //- 1/sqrt(2*Variance)
        scalar b_;
        //- Commonly used, precomputed factors
        scalar b2_, expma2b2_, abSqrtPi_, erfab_, denom_, b22denom_;
        //- Location of the peak of the PDF
        scalar m_;

    // Private member functions

        //- Newton-Raphson root finding of the CDF
        scalar newton(scalar x0, scalar y, label& nIter);

        //- Cummulative distribution function (CDF)
        inline scalar CDF(scalar x);

        //- Probability density function (PDF)
        inline scalar PDF(scalar x);

        // Disallow default bitwise copy construct and assignment
        mcInversionInletRandom(const mcInversionInletRandom&);
        void operator=(const mcInversionInletRandom&);

public:

    //- Runtime type information
        TypeName("mcInversionInletRandom");

    // Constructors

        //- Construct a random number generator given a Random object and
        //the mean velocity and the RMS of the fluctuations
        mcInversionInletRandom
        (
            Random& rnd,
            scalar Umean,
            scalar urms,
            const dictionary& dict
        );

    // Public member functions

        //- Generate one random number
        virtual scalar value();

        //- Update the coefficients for changed Umean and urms
        virtual void updateCoeffs(scalar Umean, scalar urms);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mcInversionInletRandomI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
