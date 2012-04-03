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

\*---------------------------------------------------------------------------*/

#include "mcIntegratedPositionCorrection.H"

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcIntegratedPositionCorrection, 0);
    addNamedToRunTimeSelectionTable
    (
        mcPositionCorrection,
        mcIntegratedPositionCorrection,
        mcPositionCorrection,
        integrated
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcIntegratedPositionCorrection::
mcIntegratedPositionCorrection
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& dict
)
:
    mcPositionCorrection(db, dict),

    Crho_
    (
        "Crho",
        dimless,
        lookupOrAddDefault<scalar>("Crho", 10.*mesh().time().deltaT().value())
    ),

    Cp_
    (
        "Cp",
        dimTime/dimArea,
        lookupOrAddDefault<scalar>("Cp", 10.*mesh().time().deltaT().value())
    ),

    pPosCorr_
    (
        IOobject
        (
            lookupOrDefault<word>("pPosCorrName", "pPosCorr"),
            db.time().timeName(),
            db,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),

    UPosCorr_
    (
        IOobject
        (
            lookupOrDefault<word>("UPosCorrName", "UPosCorr"),
            db.time().timeName(),
            db,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    )
{
    setRefCell
    (
        pPosCorr_,
        mesh().solutionDict().subDict("SIMPLE"),
        pRefCell_,
        pRefValue_
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::mcIntegratedPositionCorrection::correct
(
    Foam::mcParticleCloud& cloud
)
{
    const dimensionedScalar& deltaT = mesh().time().deltaT();
    const volScalarField& rho = cloud.rhocPdf();
    const volScalarField& pnd = cloud.pndcPdf();
    const volVectorField& U = cloud.Ufv();

    volScalarField rhs =
        (Crho_*rho + (1-Crho_)*pnd - 2*pnd + pnd.oldTime())/(deltaT*deltaT)
      + fvc::div(pnd*fvc::ddt(U), "posCorr::div(pmd*ddt(U))")
      + fvc::div(Crho_*U*(rho - pnd), "posCorr::div(U*(rho-pmd))")/deltaT;
    // correction factor to make the RHS of the Poisson equation integrate
    // to zero
    dimensionedScalar beta = fvc::domainIntegrate(rhs) / gSum(mesh().V());
    beta.dimensions() /= dimVolume;

    // solve Poisson equation for pPosCorr
    fvScalarMatrix pPosCorrEqn
    (
        fvm::ddt(Cp_, pPosCorr_) - fvm::laplacian(pPosCorr_) == -rhs + beta
    );
    pPosCorrEqn.solve();

    // time-integrator for correction velocity
    solve(fvm::ddt(pnd, UPosCorr_) == -fvc::grad(pPosCorr_));

    // apply
    interpolationCellPointFace<vector> UPosCorrInterp(UPosCorr_);
    forAllIter(mcParticleCloud, cloud, pIter)
    {
        mcParticle& part = pIter();
        const point& p = part.position();
        label c = part.cell();
        label f = part.face();
        part.Ucorrection() += UPosCorrInterp.interpolate(p, c, f);
    }
}

// ************************************************************************* //
