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
    Foam::mcGhostInletBoundary

Description
    Inlet boundary using ghost-layers for Monte Carlo particles

SourceFiles
    mcGhostInletBoundary.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcGhostInletBoundary_H
#define mcGhostInletBoundary_H

#include "mcGhostBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class mcGhostInletBoundary Declaration
\*---------------------------------------------------------------------------*/

class mcGhostInletBoundary :
    public mcGhostBoundary
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcGhostInletBoundary(const mcGhostInletBoundary&);
        void operator=(const mcGhostInletBoundary&);

public:

    //- Runtime type information
    TypeName("mcGhostInletBoundary");


    // Constructors

        //- Construct from components
        mcGhostInletBoundary
        (
            const fvMesh& mesh,
            label patchID,
            const dictionary& dict
        );

    // Member Functions

        //- Release ghost particles at the inlet boundary
        virtual void correct
        (
            mcParticleCloud& cloud,
            bool afterMove
        );

        //- Deal with a particle hitting the patch
        virtual void hitPatch
        (
            mcParticle& p,
            mcParticle::trackData& td
        );

        //- Deal with a particle hitting the patch without trackData
        virtual void hitPatch
        (
            mcParticle& p,
            int&
        );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
