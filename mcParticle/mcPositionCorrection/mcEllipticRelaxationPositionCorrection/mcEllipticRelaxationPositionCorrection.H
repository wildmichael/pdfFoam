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
    Foam::mcEllipticRelaxationPositionCorrection

Description
    Position correction using an elliptic relaxation model

    This is a simplified version of the algorithm described by Muradoglu et al.

    See http://dx.doi.org/10.1006/jcph.2001.6861

SourceFiles
    mcEllipticRelaxationPositionCorrection.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcEllipticRelaxationPositionCorrection_H
#define mcEllipticRelaxationPositionCorrection_H

#include "mcPositionCorrection.H"

#include "interpolationCellPointFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;

/*---------------------------------------------------------------------------*\
            Class mcEllipticRelaxationPositionCorrection Declaration
\*---------------------------------------------------------------------------*/

class mcEllipticRelaxationPositionCorrection :
    public mcPositionCorrection
{

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcEllipticRelaxationPositionCorrection
        (
            const mcEllipticRelaxationPositionCorrection&
        );
        void operator=(const mcEllipticRelaxationPositionCorrection&);

protected:

    // Protected Data

        //- Density correction potential
        volScalarField Q_;

        //- Density-correction coefficients
        dimensionedScalar kf_;
        dimensionedScalar kb_;
        dimensionedScalar N_;
        dimensionedScalar CFL_;
        dimensionedScalar eps_;
        dimensionedScalar c_;
        dimensionedScalar f_;
        dimensionedScalar b_;
        dimensionedScalar a_;

        //- Interpolator for L_
        interpolationCellPointFace<scalar> LInterp_;

public:

    //- Runtime type information
    TypeName("mcEllipticRelaxationPositionCorrection");

    // Constructors

        //- Construct from components
        mcEllipticRelaxationPositionCorrection
        (
            const objectRegistry& db,
            const dictionary& mcEllipticRelaxationPositionCorrectionDict
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
