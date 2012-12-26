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
    Foam::mcBurkeSchumannReactionModel

Description
    Burke-Schumann reaction model for infinite rate chemistry

SourceFiles
    mcBurkeSchumannReactionModel.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcBurkeSchumannReactionModel_H
#define mcBurkeSchumannReactionModel_H

#include "mcReactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class mcBurkeSchumannReactionModel Declaration
\*---------------------------------------------------------------------------*/

class mcBurkeSchumannReactionModel :
    public mcReactionModel
{

    // Private Data

        //- Specific gas constant at fuel, oxidizer and stoichiometric point
        const scalar Rfuel_, Rox_, Rstoi_;

        //- Temperature of fuel, oxidizer and stoichiometric point
        const scalar Tfuel_, Tox_, Tstoi_;

        //- Stoichiometric mixture fraction
        const scalar zstoi_;

        //- Indexes of the mixture fraction and temperature properties
        label zIdx_, TIdx_;

        //- The pressure field
        const volScalarField& p_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcBurkeSchumannReactionModel(const mcBurkeSchumannReactionModel&);
        void operator=(const mcBurkeSchumannReactionModel&);

public:

    //- Runtime type information
    TypeName("mcBurkeSchumannReactionModel");

    // Constructors

        //- Construct from components
        mcBurkeSchumannReactionModel
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
