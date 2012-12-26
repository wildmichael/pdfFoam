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
    Foam::sigStopAtWriteNow

Description
    Signal handler for interupt defined by
        OptimisationSwitches::stopAtWriteNowSignal

    Write and stop the job.

    Backported from OpenFOAM-2.1.x

SourceFiles
    sigStopAtWriteNow.C

\*---------------------------------------------------------------------------*/

#ifndef sigStopAtWriteNow_H
#define sigStopAtWriteNow_H

#ifndef FOAM_HEX_VERSION
#error FOAM_HEX_VERSION not defined
#endif

#if FOAM_HEX_VERSION < 0x200
#include <signal.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;

/*---------------------------------------------------------------------------*\
                           Class sigStopAtWriteNow Declaration
\*---------------------------------------------------------------------------*/

class sigStopAtWriteNow
{
    // Private data

        //- number of signal to use
        static int signal_;

        //- Saved old signal trapping setting
        static struct sigaction oldAction_;

    // Private Member Functions

        static void sigHandler(int);


public:

    // Constructors

        //- Construct null
        sigStopAtWriteNow();

        //- Construct from components
        sigStopAtWriteNow(const bool verbose, Time& runTime);


    //- Destructor
    ~sigStopAtWriteNow();


    // Member functions

        //- Is active?
        bool active() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // FOAM_HEX_VERSION < 0x200

#endif

// ************************************************************************* //
