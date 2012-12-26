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

Application
    initPositionCorrectionTest

Description
    Initialize the particle field for the positionCorrectionTest.

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mcParticleCloud.H"
#include "RASModel.H"
#include "basicRhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nInitializing particles\n" << endl;

    Info<< "Time = " << runTime.timeName() << nl << endl;

    cloud.clear();

    cloud.initReleaseParticles();

    {
        DimensionedField<scalar, volMesh> mMom
        (
            IOobject
            (
                "mMom",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("mMom", dimMass, 0.0)
        );

        point c(1., 1., 1.);
        label cI = mesh.findCell(c);
        forAllIter(mcParticleCloud, cloud, pIter)
        {
            mcParticle& p = pIter();
            label cellI = p.cell();
            if (cellI == cI)
            {
                p.m() *= 2.;
            }
            mMom[cellI] += p.m();
        }
        pnd.internalField() = mMom / mesh.V();
    }

    cloud.write();
    pnd.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
