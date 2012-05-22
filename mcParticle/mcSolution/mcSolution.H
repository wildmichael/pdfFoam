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
    Foam::mcSolution

Description
    Selector class for relaxation factors, models and their coefficients.

SourceFiles
    mcSolution.C

\*---------------------------------------------------------------------------*/

#ifndef mcSolution_H
#define mcSolution_H

#include "IOdictionary.H"
#include "dimensionedScalar.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class mcSolution Declaration
\*---------------------------------------------------------------------------*/

class mcSolution:
    public IOdictionary
{
    // Private data

        //- The averaging time scale
        dimensionedScalar averagingTime_;
        //- Default relaxation time scale
        dimensionedScalar defaultRelaxationTime_;
        //- The relaxation time scales
        dictionary relaxationTimes_;
        //- The number of particles per cell
        label particlesPerCell_;
        //- Whether particle number control is enabled
        Switch particleNumberControl_;
        //- Cloning threshold
        scalar cloneAt_;
        //- Elimination threshold
        scalar eliminateAt_;
        //- Minimum TKE
        dimensionedScalar kMin_;
        //- Numerical diffusion coefficient
        dimensionedScalar DNum_;

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        mcSolution(const mcSolution&);
        void operator=(const mcSolution&);


public:


    // Constructors

        //- Construct for given objectRegistry and dictionary
        mcSolution(const objectRegistry& obr);


    // Member Functions

        // Access

            //- Return the selected sub-dictionary if the "select" keyword is
            //  given, otherwise return the complete dictionary
            const dictionary& solutionDict() const;

            //- Return the averaging time scale
            const dimensionedScalar& averagingTime() const
            {return averagingTime_;}

            //- Return the relaxation time for the given field
            dimensionedScalar relaxationTime(const word& name) const;

            //- Return the number of particles per cell
            label particlesPerCell() const {return particlesPerCell_;}

            //- Return whether particle number control is enabled
            bool enableParticleNumberControl() const
            {return particleNumberControl_;}

            //- Return the cloning threshold
            const scalar& cloneAt() const {return cloneAt_;}

            //- Return the elimination threshold
            const scalar& eliminateAt() const {return eliminateAt_;}

            //- Return the minimum TKE
            const dimensionedScalar& kMin() const {return kMin_;}

            //- Return the numerical diffusion coefficient
            const dimensionedScalar& DNum() const {return DNum_;}

        // Read

            //- Read the mcSolution dictionary
            bool read();

            //- Read object if modified
            bool readIfModified();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
