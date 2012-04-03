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
    Foam::mcIntegratedPositionCorrection

Description
    Position correction scheme using time-integration

    Reasoning:
    - continuity:
        ddt(pmd) + div((Umean + Uc)*pmd) == 0                               (1)
    - time integrator:
        ddt(pmd, Uc) == -grad(pc)                                           (2)
    - take divergence of (2), substitute (1):
        laplacian(pc) == d2dt2(pmd) + div(ddt(pmd, Umean))
                      == d2dt2(pmd) + div(pmd*ddt(Umean) + Umean*ddt(pmd))  (3)
    - assuming that pmd(t+deltaT) == rho:
        laplacian(pc) ==  (rho - 2*pmd + pmd.oldTime())/(deltaT*deltaT)
                        + div(pmd*ddt(Umean) + Umean*(rho - pmd)/deltaT)    (4)
    - underrelax rho as (Crho*rho + (1-Crho)*pmd), add compressibility Cp
      to (4):
        ddt(Cp,pc) == laplacian(pc)
                    - (Crho*rho - (1+Crho)*pmd + pmd.oldtime())/(deltaT*deltaT)
                    + div(pmd*ddt(Umean) + Crho*Umean*(rho - pmd)/deltaT)   (5)
    - solve (5), then (2)

SourceFiles
    mcIntegratedPositionCorrection.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcIntegratedPositionCorrection_H
#define mcIntegratedPositionCorrection_H

#include "mcPositionCorrection.H"

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;

/*---------------------------------------------------------------------------*\
                Class mcIntegratedPositionCorrection Declaration
\*---------------------------------------------------------------------------*/

class mcIntegratedPositionCorrection :
    public mcPositionCorrection
{
    // Private Data

        //- Density-correction coefficients
        dimensionedScalar Crho_;
        dimensionedScalar Cp_;

        //- Correction pressure
        volScalarField pPosCorr_;

        //- Reference cell
        label pRefCell_;

        //- Reference pressure
        scalar pRefValue_;

        //- Correction velocity
        volVectorField UPosCorr_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcIntegratedPositionCorrection(const mcIntegratedPositionCorrection&);
        void operator=(const mcIntegratedPositionCorrection&);

public:

    //- Runtime type information
    TypeName("mcIntegratedPositionCorrection");

    // Constructors

        //- Construct from components
        mcIntegratedPositionCorrection
        (
            const objectRegistry& db,
            const dictionary& mcIntegratedPositionCorrectionDict
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
