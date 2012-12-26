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
    Foam::mcIEMMixingModel

Description
    Interaction by exchange with the mean (IEM) mixing model

SourceFiles
    mcIEMMixingModel.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcIEMMixingModel_H
#define mcIEMMixingModel_H

#include "autoPtr.H"
#include "dictionary.H"
#include "mcMixingModel.H"
#include "runTimeSelectionTables.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;

/*---------------------------------------------------------------------------*\
                        Class mcIEMMixingModel Declaration
\*---------------------------------------------------------------------------*/

class mcIEMMixingModel :
    public mcMixingModel
{

    // Private Data

        //- IEM mixing model constant Cmix / 2
        scalar Cmix2_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcIEMMixingModel(const mcIEMMixingModel&);
        void operator=(const mcIEMMixingModel&);

public:

    //- Runtime type information
    TypeName("mcIEMMixingModel");


    // Constructors

        //- Construct from components
        mcIEMMixingModel
        (
            mcParticleCloud& cloud,
            const objectRegistry& db,
            const word& subDictName
        );


    // Member Functions

        //- Update/prepare internal data
        virtual void updateInternals();

        //- Apply the mixing model to a single particle
        virtual void correct(mcParticle& p);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
