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
    Foam::mcGhostInletOutletBoundary

Description
    Inlet/outlet boundary using ghost-layers for Monte Carlo particles

SourceFiles
    mcGhostInletOutletBoundary.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcGhostInletOutletBoundary_H
#define mcGhostInletOutletBoundary_H

#include "mcOpenBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class mcGhostInletOutletBoundary Declaration
\*---------------------------------------------------------------------------*/

class mcGhostInletOutletBoundary :
    public mcOpenBoundary
{
    // Private data

        //- Global hash of ghost cell indices
        static HashSet<label> ghostCellHash_;
        //- Global indicator whether ghost-particles have been purged
        static bool purgedGhosts_;
        //- List of cell-indices for each ghost patch
        labelList    ghostCellLayer_;
        //- List of ghost-cell shifts for each ghost patch
        vectorList ghostCellShifts_;
#ifdef FULLDEBUG
        //- Output for generated particle velocities (only those
        // <em>entering</em> the domain)
        autoPtr<OFstream> statFile_;
#endif

    // Private Member Functions

        //- Enforce in/out flow BCs by populating ghost cells
        void populateGhostCells();

        //- Remove all ghost particles from cloud
        void purgeGhostParticles();

        //- Finds the ghost cell layers
        void findGhostLayer();

        // Disallow default bitwise copy construct and assignment
        mcGhostInletOutletBoundary(const mcGhostInletOutletBoundary&);
        void operator=(const mcGhostInletOutletBoundary&);

public:

    //- Runtime type information
    TypeName("mcGhostInletOutletBoundary");


    // Constructors

        //- Construct from components
        mcGhostInletOutletBoundary
        (
            mcParticleCloud& cloud,
            label patchID,
            const dictionary& dict
        );

    // Member Functions

        //- Release ghost particles at the boundary
        virtual void correct(bool afterMove);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
