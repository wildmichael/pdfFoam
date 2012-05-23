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
    Foam::mcNaiveFlameletModel

Description
    Dummy reaction model for diffusion flame simulations

SourceFiles
    mcNaiveFlameletModel.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcNaiveFlameletModel_H
#define mcNaiveFlameletModel_H

#include "mcReactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class mcNaiveFlameletModel Declaration
\*---------------------------------------------------------------------------*/

class mcNaiveFlameletModel :
    public mcReactionModel
{

    // Private Data

        label zIdx_, TIdx_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcNaiveFlameletModel(const mcNaiveFlameletModel&);
        void operator=(const mcNaiveFlameletModel&);

public:

    //- Runtime type information
    TypeName("mcNaiveFlameletModel");

    // Constructors

        //- Construct from components
        mcNaiveFlameletModel
        (
            mcParticleCloud& cloud,
            const objectRegistry& db,
            const word& subDictName
        );

    // Member Functions

        //- Apply the reaction model to a single particle
        virtual void correct(mcParticle& p);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
