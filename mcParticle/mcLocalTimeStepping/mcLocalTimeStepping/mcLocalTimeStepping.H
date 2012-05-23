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
    Foam::mcLocalTimeStepping

Description
    Abstract base class for local time-stepping models

SourceFiles
    mcLocalTimeStepping.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcLocalTimeStepping_H
#define mcLocalTimeStepping_H

#include "mcModel.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;
class mcParticle;

/*---------------------------------------------------------------------------*\
                     Class mcLocalTimeStepping Declaration
\*---------------------------------------------------------------------------*/

class mcLocalTimeStepping :
    public mcModel
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcLocalTimeStepping(const mcLocalTimeStepping&);
        void operator=(const mcLocalTimeStepping&);

public:

    //- Runtime type information
    TypeName("mcLocalTimeStepping");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mcLocalTimeStepping,
            mcLocalTimeStepping,
            (
                mcParticleCloud& cloud,
                const objectRegistry& db,
                const word& subDictName
            ),
            (cloud, db, subDictName)
        );

    // Constructors

        //- Construct from components
        mcLocalTimeStepping
        (
            mcParticleCloud& cloud,
            const objectRegistry& db,
            const word& subDictName
        );


    // Selectors

        //- Return a reference to the selected position correction
        static autoPtr<mcLocalTimeStepping> New
        (
            mcParticleCloud& cloud,
            const objectRegistry& db
        );

    // Member Functions

        using mcModel::correct;

        //- Apply the position correction to a single particle
        virtual void correct(mcParticle& p);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
