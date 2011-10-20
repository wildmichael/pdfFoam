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
    Foam::mcInletOutletBoundary

Description
    Inlet/outlet boundary for Monte Carlo particles

SourceFiles
    mcInletOutletBoundary.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcInletOutletBoundary_H
#define mcInletOutletBoundary_H

#include "mcOpenBoundary.H"

#include "InletRandom.H"

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

        //- RMS of fluctuations at boundary in face-normal coordinate system
        vectorField urms_;

        //- Random number generators for the wall-normal velocity
        // Do not access directly, use getInRnd() instead.
        PtrList<InletRandom> inRnd_;

#ifdef FULLDEBUG
        //- Output for generated particle velocities
        autoPtr<OFstream> statFile_;
#endif

    // Private Member Functions

        //- Get the U_ data field in face-normal coordinate system and create
        // it if necessary
        const vectorField& getU(const mcParticleCloud& cloud);

        //- Get the inRnd_ list and create it if necessary
        PtrList<InletRandom>& getInRnd(mcParticleCloud& cloud);

        //- Generate random point on a face
        point randomPoint(mcParticleCloud& cloud, label faceI);

        //- Generate random velocity on a face
        vector randomVelocity(mcParticleCloud& cloud, label faceI);

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
            const fvMesh& mesh,
            label patchID,
            const dictionary& dict
        );

    // Member Functions

        //- Release particles at the inlet boundary
        virtual void correct
        (
            mcParticleCloud& cloud,
            bool afterMove
        );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
