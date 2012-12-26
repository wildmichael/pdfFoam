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

        //- Location of the peak of the PDF
        scalar m_;

    // Private member functions

        //- Newton-Raphson root finding of the CDF
        scalar newton(scalar x0, scalar y, label& nIter);

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
            scalar UMean,
            scalar uRms,
            const dictionary& dict
        );

    // Public member functions

        //- Generate one random number
        virtual scalar value();

        //- Update the coefficients for changed UMean and uRms
        virtual void updateCoeffs(scalar UMean, scalar uRms);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
