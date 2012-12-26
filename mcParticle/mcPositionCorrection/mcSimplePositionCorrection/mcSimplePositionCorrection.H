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
    Foam::mcSimplePositionCorrection

Description
    Simple position correction model

    Uses the face area of the cell bounding box and the magnitude of the mean
    velocity to scale the correction velocity. The idea behind the former is
    that the mass exchange should be independent of the face area. The reasoning
    for the latter is that the transferred momentum during deltaT should be
    independent of the transport.

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
template<class> class interpolation;

/*---------------------------------------------------------------------------*\
                  Class mcSimplePositionCorrection Declaration
\*---------------------------------------------------------------------------*/

class mcSimplePositionCorrection :
    public mcPositionCorrection
{
    // Private Data

        //- Cell-bounding-box cross-sectional inverse areas
        DimensionedField<vector, volMesh> Ainv_;

        //- Correction potential
        volScalarField phi_;

        //- Gradient of the correction potential
        volVectorField gradPhi_;

        //- Interpolator for the correction potential
        autoPtr<interpolation<vector> > gradPhiInterp_;

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
