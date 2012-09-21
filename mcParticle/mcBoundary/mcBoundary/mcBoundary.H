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
    Foam::mcBoundary

Description
    Abstract base class for Monte Carlo boundary condition handlers

SourceFiles
    mcBoundaryI.H
    mcBoundary.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcBoundary_H
#define mcBoundary_H

#include "dictionary.H"
#include "typeInfo.H"
#include "fvMesh.H"
#include "mcParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;

/*---------------------------------------------------------------------------*\
                          Class mcBoundary Declaration
\*---------------------------------------------------------------------------*/

class mcBoundary :
    public dictionary
{

    // Private Data

        //- The particle cloud
        mcParticleCloud& cloud_;

        //- The associated patch ID
        const label patchID_;

        //- The associated patch
        const polyPatch& patch_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcBoundary(const mcBoundary&);
        void operator=(const mcBoundary&);

public:

    //- Runtime type information
    TypeName("mcBoundary");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mcBoundary,
            mcBoundary,
            (
                mcParticleCloud& cloud,
                label patchID,
                const dictionary& dict
            ),
            (cloud, patchID, dict)
        );


    // Constructors

        //- Construct from components
        mcBoundary
        (
            mcParticleCloud& cloud,
            label patchID,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected reaction model
        static autoPtr<mcBoundary> New
        (
            mcParticleCloud& cloud,
            label patchID,
            const dictionary& dict
        );

    //- Destructor
        virtual ~mcBoundary()
        {}

    // Member Functions

        //- Apply the boundary condition to the particle cloud, does nothing by
        // default
        virtual void correct(bool afterMove)
        {}

        //- Deal with a particle hitting the patch
        virtual
        void hitPatch
        (
            mcParticle& p,
            mcParticle::trackData& td
        ) = 0;

        //- Deal with a particle hitting the patch withou trackData
        virtual
        void hitPatch
        (
            mcParticle& p,
            int&
        ) = 0;

        //- The particle cloud
        inline const mcParticleCloud& cloud() const;

        //- The particle cloud
        inline mcParticleCloud& cloud();

        //- The ID of the associated patch
        inline label patchID() const;

        inline const polyPatch& patch() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mcBoundaryI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
