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
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;
template<class> class interpolation;

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
        //- Normalized density difference
        volScalarField QInst_;
        //- Gradients of the density difference and correction potential
        volVectorField gradQInst_, gradQ_;
        //- Indicator variable
        volScalarField zeta_;
        // maximum mean velocity
        dimensionedScalar U0_;

        //- Density-correction coefficients
        dimensionedScalar eps_;
        dimensionedScalar c_;
        dimensionedScalar f_;
        dimensionedScalar b_;
        dimensionedScalar a_;

        //- Interpolators
        autoPtr<interpolation<scalar> > LInterp_;
        // TODO try gradInterpolationConstantTet
        autoPtr<interpolation<vector> > gradQInterp_;
        autoPtr<interpolation<vector> > gradQInstInterp_;
        autoPtr<interpolation<scalar> > zetaInterp_;

    // Protected Member Functions

        //- Update coefficients
        void readCoeffs();

public:

    //- Runtime type information
    TypeName("mcEllipticRelaxationPositionCorrection");

    // Constructors

        //- Construct from components
        mcEllipticRelaxationPositionCorrection
        (
            mcParticleCloud& cloud,
            const objectRegistry& db,
            const word& subDictName
        );


    // Member Functions

        //- Update/prepare internal data
        virtual void updateInternals();

        //- Apply the position correction
        virtual void correct(mcParticle& p);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
