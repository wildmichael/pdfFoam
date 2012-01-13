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
    Foam::mcParticleLocalTimeStepping

Description
    Cell-based local time-stepping

    As described in Hegetschweiler, 2010 (doi:10.3929/ethz-a-006409982).
    Note that there are various error in mentioned document.

SourceFiles
    mcParticleLocalTimeStepping.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcParticleLocalTimeStepping_H
#define mcParticleLocalTimeStepping_H

#include "mcLocalTimeStepping.H"

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class mcParticleLocalTimeStepping Declaration
\*---------------------------------------------------------------------------*/

class mcParticleLocalTimeStepping :
    public mcLocalTimeStepping
{
    // Private Data

        //- Courant number
        scalar CourantU_;

    // Private Member Functions

        //- Update the eta field
        void updateEta(const mcParticleCloud& cloud);

        // Disallow default bitwise copy construct and assignment
        mcParticleLocalTimeStepping(const mcParticleLocalTimeStepping&);
        void operator=(const mcParticleLocalTimeStepping&);

public:

    //- Runtime type information
    TypeName("mcParticleLocalTimeStepping");

    // Constructors

        //- Construct from components
        mcParticleLocalTimeStepping
        (
            const Foam::objectRegistry& db,
            const dictionary& dict
        );

    // Member Functions

        //- Apply the position correction
        virtual void correct(mcParticleCloud& cloud);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
