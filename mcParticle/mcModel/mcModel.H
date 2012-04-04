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
    Foam::mcModel

Description
    Abstract base class for Monte Carlo sub-models

SourceFiles
    mcModelI.H
    mcModel.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcModel_H
#define mcModel_H

#include "dictionary.H"
#include "labelList.H"
#include "List.H"
#include "typeInfo.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;
class objectRegistry;

/*---------------------------------------------------------------------------*\
                           Class mcModel Declaration
\*---------------------------------------------------------------------------*/

class mcModel :
    public dictionary
{

    // Private Data

        //- Object registry
        const objectRegistry& db_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcModel(const mcModel&);
        void operator=(const mcModel&);

public:

    //- Runtime type information
    TypeName("mcModel");

    // Constructors

        //- Construct from components
        mcModel
        (
            const objectRegistry& db,
            const dictionary& dict
        );

    //- Destructor
        virtual ~mcModel()
        {}

    // Member Functions

        //- Apply the model to the particle cloud
        virtual void correct(mcParticleCloud& cloud) = 0;

        //- The object registry
        inline const objectRegistry& db() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mcModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
