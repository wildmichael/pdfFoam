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
    Foam::mcInletOutletBoundary

Description
    Inlet/outlet boundary for Monte Carlo particles

    Mean velocity (getU()) and Reynolds Averaged Stress (getR()) are
    initialized from the U and TauCloudPDF fields, respectively. However,
    currently only the isotropic part of TauCloudPDF is used. Generating
    correlated fluctuating velocities is on the to-do list...

ToDo
    Use the deviatoric part of TauCloudPDF to determine the velocity
    fluctuation correlations and generate random fluctuations accordingly.

SourceFiles
    mcInletOutletBoundary.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcInletOutletBoundary_H
#define mcInletOutletBoundary_H

#include "mcOpenBoundary.H"

#include "mcInletRandom.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class mcInletOutletBoundary Declaration
\*---------------------------------------------------------------------------*/

class Random;

class mcInletOutletBoundary :
    public mcOpenBoundary
{

    // Private data

        //- Weights of surface triangles of all faces
        List<scalarList> probVec_;

        //- Transformations from global into face-normal system
        tensorField fwdTrans_;

        //- Transformations from face-normal into global system
        tensorField revTrans_;

        //- Mean velocity in face-normal coordinate system
        // Do not access directly, use getU() instead.
        vectorField U_;

        //- RAS at boundary in face-normal coordinate system
        // Do not access directly, use getR() instead.
        symmTensorField R_;

        //- Random number generators for the wall-normal velocity
        // Do not access directly, use getInRnd() instead.
        PtrList<mcInletRandom> inRnd_;

#ifdef FULLDEBUG
        //- Output for generated particle velocities
        autoPtr<OFstream> statFile_;
#endif

    // Private Member Functions

        //- Get the U_ data field in face-normal coordinate system and create
        // it if necessary
        const vectorField& getU();

        //- Get the R_ data field in face-normal coordinate system and create
        // it if necessary
        const symmTensorField& getR();

        //- Get the inRnd_ list and create it if necessary
        PtrList<mcInletRandom>& getInRnd();

        //- Generate random point on a face
        point randomPoint(label faceI);

        //- Generate random velocity on a face
        vector randomVelocity(label faceI);

        // Disallow default bitwise copy construct and assignment
        mcInletOutletBoundary(const mcInletOutletBoundary&);
        void operator=(const mcInletOutletBoundary&);

public:

    //- Runtime type information
    TypeName("mcInletOutletBoundary");


    // Constructors

        //- Construct from components
        mcInletOutletBoundary
        (
            mcParticleCloud& cloud,
            label patchID,
            const dictionary& dict
        );

    // Member Functions

        //- Release particles at the inlet boundary
        virtual void correct(bool afterMove);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
