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
    Foam::mcGhostBoundary

Description
    Base-class for Monte Carlo particle boundaries that require ghost layers

SourceFiles
    mcGhostBoundary.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcGhostBoundary_H
#define mcGhostBoundary_H

#include "mcBoundary.H"

#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mcGhostBoundary Declaration
\*---------------------------------------------------------------------------*/

class mcGhostBoundary :
    public mcBoundary
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

    // Private Member Functions

        //- Enforce in/out flow BCs by populating ghost cells
        void populateGhostCells(mcParticleCloud& cloud);

        //- Remove all ghost particles from cloud
        void purgeGhostParticles(mcParticleCloud& cloud);

        //- Finds the ghost cell layers
        void findGhostLayer();

        // Disallow default bitwise copy construct and assignment
        mcGhostBoundary(const mcGhostBoundary&);
        void operator=(const mcGhostBoundary&);

public:

    //- Runtime type information
    TypeName("mcGhostBoundary");


    // Constructors

        //- Construct from components
        mcGhostBoundary
        (
            const fvMesh& mesh,
            label patchID,
            const dictionary& dict
        );

    // Member Functions

        //- Release particles at the inlet boundary
        virtual void correct(mcParticleCloud& cloud, bool afterMove);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
