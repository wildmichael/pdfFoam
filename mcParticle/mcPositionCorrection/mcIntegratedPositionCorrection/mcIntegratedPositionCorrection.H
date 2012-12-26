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
    Foam::mcIntegratedPositionCorrection

Description
    Position correction scheme using a global Poisson equation to find a
    correction pressure whose negative gradient is integrated in time to obtain
    the correction velocity.

SourceFiles
    mcIntegratedPositionCorrection.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcIntegratedPositionCorrection_H
#define mcIntegratedPositionCorrection_H

#include "mcPositionCorrection.H"

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;
template<class> class interpolation;

/*---------------------------------------------------------------------------*\
                Class mcIntegratedPositionCorrection Declaration
\*---------------------------------------------------------------------------*/

class mcIntegratedPositionCorrection :
    public mcPositionCorrection
{
    // Private Data

        //- Correction pressure
        volScalarField pPosCorr_;

        //- Reference cell
        label pRefCell_;

        //- Reference pressure
        scalar pRefValue_;

        //- Correction velocity
        volVectorField UPosCorr_;

        //- Correction velocity interpolator
        autoPtr<interpolation<vector> > UPosCorrInterp_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcIntegratedPositionCorrection(const mcIntegratedPositionCorrection&);
        void operator=(const mcIntegratedPositionCorrection&);

public:

    //- Runtime type information
    TypeName("mcIntegratedPositionCorrection");

    // Constructors

        //- Construct from components
        mcIntegratedPositionCorrection
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
