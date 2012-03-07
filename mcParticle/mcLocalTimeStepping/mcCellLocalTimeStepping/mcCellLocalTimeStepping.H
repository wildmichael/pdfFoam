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

#include "interpolationCellPointFace.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class mcCellLocalTimeStepping Declaration
\*---------------------------------------------------------------------------*/

class mcCellLocalTimeStepping :
    public mcLocalTimeStepping
{
    // Private Data

        //- Courant number
        scalar CourantU_;

        //- Maximum local to global time step ratio
        scalar etaMax_;

        //- Local to global time step ratio field
        volScalarField eta_;

        //- Interpolator for Omega
        autoPtr<interpolationCellPointFace<scalar> > etaInterp_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcCellLocalTimeStepping(const mcCellLocalTimeStepping&);
        void operator=(const mcCellLocalTimeStepping&);

protected:

    // Protected Member Functions

        //- Update the eta field and the interpolator
        virtual void setupInternals(const mcParticleCloud& cloud);

public:

    //- Runtime type information
    TypeName("mcCellLocalTimeStepping");

    // Constructors

        //- Construct from components
        mcCellLocalTimeStepping
        (
            const Foam::objectRegistry& db,
            const dictionary& dict
        );

    // Member Functions

        //- Apply the position correction to a single particle
        virtual void correct
        (
            mcParticleCloud& cloud,
            mcParticle& p,
            bool prepare=true
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
