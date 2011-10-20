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
    Foam::mcInletRandom

Description
    Abstract base class for inlet velocity random number generators

    Generate random values from the PDF
    \f[ f_x = \frac{X}{N} e^{-(X-a)^2 b^2} \f]
    where
    \f{eqnarray*}{
        b &=& \frac{1}{\sqrt{2}u_\mathrm{rms}} \quad \text{and} \\
        N &=& b^2 e^{-\frac{a^2}{2 b^2}} + a \sqrt{\frac{\pi b^2}{2}}
              \left(1+\mathrm{erf}\left(\frac{a}{\sqrt{2 b^2}}\right)\right)
    \f}
    using various methods.


SourceFiles
    mcInletRandom.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcInletRandom_H
#define mcInletRandom_H

#include "dictionary.H"
#include "scalar.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Random;

/*---------------------------------------------------------------------------*\
                        Class mcInletRandom Declaration
\*---------------------------------------------------------------------------*/

class mcInletRandom
    : public dictionary
{
    // Private data

        //- Random number generator
        Random& rnd_;
        //- Mean velocity
        scalar Umean_;
        //- RMS of fluctuations (stddev)
        scalar urms_;

    // Private member functions

        // Disallow default bitwise copy construct and assignment
        mcInletRandom(const mcInletRandom&);
        void operator=(const mcInletRandom&);

public:

    //- Runtime type information
    TypeName("mcInletRandom");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mcInletRandom,
            mcInletRandom,
            (
                Random& rnd,
                scalar Umean,
                scalar urms,
                const dictionary& dict
            ),
            (rnd, Umean, urms, dict)
        );

    // Constructors

        //- Construct a random number generator given a Random object and
        // the mean velocity and the RMS of the fluctuations
        mcInletRandom
        (
            Random& rnd,
            scalar Umean,
            scalar urms,
            const dictionary& dict
        );

    //- Destructor

        virtual ~mcInletRandom();


    // Selectors

        //- Return a reference to the selected random number generator
        static autoPtr<mcInletRandom> New
        (
            Random& rnd,
            scalar Umean,
            scalar urms,
            const dictionary& dict
        );

    // Access

        Random& rnd() { return rnd_; }
        scalar Umean() const { return Umean_; }
        scalar& Umean() { return Umean_; }
        scalar urms() const { return urms_; }
        scalar& urms() { return urms_; }

    // Public member functions

        //- Generate one random number
        virtual scalar value() = 0;

        //- Update the coefficients for changed Umean and urms
        virtual void updateCoeffs(scalar Umean, scalar urms);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
