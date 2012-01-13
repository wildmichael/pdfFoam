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
    Foam::mcSimplePositionCorrection

Description
    Simple position correction model

    Uses the face area of the cell bounding box and the magnitude of the mean
    velocity to scale the correction velocity.

SourceFiles
    mcSimplePositionCorrection.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcSimplePositionCorrection_H
#define mcSimplePositionCorrection_H

#include "mcPositionCorrection.H"

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;

/*---------------------------------------------------------------------------*\
                  Class mcSimplePositionCorrection Declaration
\*---------------------------------------------------------------------------*/

class mcSimplePositionCorrection :
    public mcPositionCorrection
{
    // Private Data

        //- Density-correction coefficients
        dimensionedScalar C_;

        //- Cell-bounding-box cross-sectional inverse areas
        DimensionedField<vector, volMesh> Ainv_;

        //- Correction potential
        volScalarField phi_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcSimplePositionCorrection(const mcSimplePositionCorrection&);
        void operator=(const mcSimplePositionCorrection&);

public:

    //- Runtime type information
    TypeName("mcSimplePositionCorrection");

    // Constructors

        //- Construct from components
        mcSimplePositionCorrection
        (
            const objectRegistry& db,
            const dictionary& mcSimplePositionCorrectionDict
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
