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
    Foam::mcVelocityModel

Description
    Abstract base class for Monte Carlo particle velocity models

SourceFiles
    mcVelocityModel.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcVelocityModel_H
#define mcVelocityModel_H

#include "mcModel.H"

#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "mcParticle.H"
#include "runTimeSelectionTables.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mcVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class mcVelocityModel :
    public mcModel
{

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcVelocityModel(const mcVelocityModel&);
        void operator=(const mcVelocityModel&);

public:

    //- Runtime type information
    TypeName("mcVelocityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mcVelocityModel,
            mcVelocityModel,
            (
                mcParticleCloud& cloud,
                const objectRegistry& db,
                const word& subDictName
            ),
            (cloud, db, subDictName)
        );

    // Constructors

        //- Construct from components
        mcVelocityModel
        (
            mcParticleCloud& cloud,
            const objectRegistry& db,
            const word& subDictName
        );


    // Selectors

        //- Return a reference to the selected velocity model
        static autoPtr<mcVelocityModel> New
        (
            mcParticleCloud& cloud,
            const objectRegistry& db
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
