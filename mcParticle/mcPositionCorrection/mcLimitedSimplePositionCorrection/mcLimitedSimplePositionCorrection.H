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
    Foam::mcLimitedSimplePositionCorrection

Description
    Limited, simple position correction model

    As used in PDF2DFV

SourceFiles
    mcLimitedSimplePositionCorrection.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcLimitedSimplePositionCorrection_H
#define mcLimitedSimplePositionCorrection_H

#include "mcPositionCorrection.H"

#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;

/*---------------------------------------------------------------------------*\
              Class mcLimitedSimplePositionCorrection Declaration
\*---------------------------------------------------------------------------*/

class mcLimitedSimplePositionCorrection :
    public mcPositionCorrection
{
    // Private Data

        //- Density-correction coefficients
        dimensionedScalar C_;
        dimensionedScalar CFL_;
        dimensionedScalar Cd_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcLimitedSimplePositionCorrection
        (
            const mcLimitedSimplePositionCorrection&
        );
        void operator=(const mcLimitedSimplePositionCorrection&);

public:

    //- Runtime type information
    TypeName("mcLimitedSimplePositionCorrection");

    // Constructors

        //- Construct from components
        mcLimitedSimplePositionCorrection
        (
            const objectRegistry& db,
            const dictionary& parentDict,
            const dictionary& mcLimitedSimplePositionCorrectionDict
        );


    // Member Functions

        //- Apply the position correction
        virtual void correct(mcParticleCloud& cloud);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
