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

class mcParticle;
class mcParticleCloud;
class objectRegistry;

/*---------------------------------------------------------------------------*\
                           Class mcModel Declaration
\*---------------------------------------------------------------------------*/

class mcModel
{

    // Private Data

        //- Particle cloud
        mcParticleCloud& cloud_;

        //- Object registry
        const objectRegistry& db_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcModel(const mcModel&);
        void operator=(const mcModel&);

protected:

    // Protected Member Data

        //- Name of the constant/mcSolution subdict
        word subDictName_;

        //- Dictionary in constant/thermophysicalProperties
        dictionary thermoDict_;

public:

    //- Runtime type information
    TypeName("mcModel");

    // Constructors

        //- Construct from components
        mcModel
        (
            mcParticleCloud& cloud,
            const objectRegistry& db,
            const word& subDictName
        );

    //- Destructor
        virtual ~mcModel()
        {}

    // Member Functions

        //- The constant/thermophysicalModels subdict
        inline const dictionary& thermoDict() const;

        //- The system/mcSolution subdict
        // @note If the model-specific sub-dictionary does not exist in the
        // sub-dictionary specified in a @c select entry, the top-level
        // sub-dictionary will be returned. This is to avoid reverting back
        // to the default settings.
        dictionary solutionDict() const;

        //- Update/prepare internal data
        virtual void updateInternals();

        //- Apply the model to the particle cloud
        virtual void correct();

        //- Apply the model to a single particle
        // @note It is the callers responsibility to call updateInternals
        // if required before calling this function.
        virtual void correct(mcParticle& p) = 0;

        //- Compute the Courant number due to this model for given particle.
        // Unless overridden does nothing.
        virtual void Co(mcParticle&) const;

        //- The particle cloud
        inline const mcParticleCloud& cloud() const;

        //- The particle cloud
        inline mcParticleCloud& cloud();

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
