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
    Foam::mcEmptyBoundary

Description
    Empty boundary for Monte Carlo particles. It is a fatal error if a particle
    hits this boundary.

SourceFiles
    mcEmptyBoundary.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcEmptyBoundary_H
#define mcEmptyBoundary_H

#include "mcBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mcEmptyBoundary Declaration
\*---------------------------------------------------------------------------*/

class mcEmptyBoundary :
    public mcBoundary
{

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcEmptyBoundary(const mcEmptyBoundary&);
        void operator=(const mcEmptyBoundary&);

public:

    //- Runtime type information
    TypeName("mcEmptyBoundary");


    // Constructors

        //- Construct from components
        mcEmptyBoundary
        (
            mcParticleCloud& cloud,
            label patchID,
            const dictionary& dict
        );

    // Member Functions

        //- Deal with a particle hitting the patch
        virtual void hitPatch
        (
            mcParticle& p,
#if FOAM_HEX_VERSION < 0x200
            mcParticle::trackData& td
#else
            mcParticle::trackData& td,
            const label patchI,
            const scalar trackFraction,
            const tetIndices& tetIs
#endif
        );

#if FOAM_HEX_VERSION < 0x200
        //- Deal with a particle hitting the patch without trackData
        virtual void hitPatch
        (
            mcParticle& p,
            int&
        );
#endif

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
