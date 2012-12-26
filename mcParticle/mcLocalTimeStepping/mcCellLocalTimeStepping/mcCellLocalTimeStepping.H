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
    Foam::mcCellLocalTimeStepping

Description
    Cell-based local time-stepping

    As described in Hegetschweiler, 2010 (doi:10.3929/ethz-a-006409982).
    Note that there are various error in mentioned document.

SourceFiles
    mcCellLocalTimeStepping.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcCellLocalTimeStepping_H
#define mcCellLocalTimeStepping_H

#include "mcLocalTimeStepping.H"

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class> class interpolation;

/*---------------------------------------------------------------------------*\
                   Class mcCellLocalTimeStepping Declaration
\*---------------------------------------------------------------------------*/

class mcCellLocalTimeStepping :
    public mcLocalTimeStepping
{
    // Private Data

        //- Local to global time step ratio field
        volScalarField eta_;

        //- Interpolator for Omega
        autoPtr<interpolation<scalar> > etaInterp_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcCellLocalTimeStepping(const mcCellLocalTimeStepping&);
        void operator=(const mcCellLocalTimeStepping&);

public:

    //- Runtime type information
    TypeName("mcCellLocalTimeStepping");

    // Constructors

        //- Construct from components
        mcCellLocalTimeStepping
        (
            mcParticleCloud& cloud,
            const objectRegistry& db,
            const word& subDictName
        );

    // Member Functions

        //- Update the eta field and the interpolator
        virtual void updateInternals();

        //- Apply the position correction to a single particle
        virtual void correct(mcParticle& p);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
