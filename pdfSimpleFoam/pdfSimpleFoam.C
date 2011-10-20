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

Application
    rhoSimpleFoam

Description
    Steady-state SIMPLE solver for laminar or turbulent RANS flow of
    compressible fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mcThermo.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    bool prevCycleWasFV = false;
    volVectorField gradP = fvc::grad(p);

    scalar eqnResidual = 1, maxFVResidual = 0, maxPDFResidual = 0;
    scalar FVConvergenceCriterion = 0, PDFConvergenceCriterion = 0;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readThermoControls.H"
        #include "readSIMPLEControls.H"

        simple.readIfPresent("convergence", FVConvergenceCriterion);
        thermo.readIfPresent("convergence", PDFConvergenceCriterion);

        if (FVCycle)
        {
            maxFVResidual = 0.;

            p.storePrevIter();
            rho.storePrevIter();

            // Pressure-velocity SIMPLE corrector
            {
                // TODO bound rho?
                rho.relax();
                Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;

                #include "UEqn.H"
                #include "pEqn.H"
            }
            turbulence->correct();
            prevCycleWasFV = true;
        }
        else
        {
            maxPDFResidual = 0.;

            if (prevCycleWasFV)
            {
                gradP = fvc::grad(p);
            }
            eqnResidual = thermo.evolve();
            maxPDFResidual = max(eqnResidual, maxPDFResidual);
            prevCycleWasFV = false;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        #include "convergenceCheck.H"
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
