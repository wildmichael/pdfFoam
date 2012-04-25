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
    const Foam::dictionary& parentDict,
    const Foam::dictionary& mcIntegratedPositionCorrectionDict
)
:
    mcPositionCorrection(db, parentDict, mcIntegratedPositionCorrectionDict),

    C_
    (
        "C",
        dimless,
        lookupOrAddDefault<scalar>("C", 1e-5)
    ),

    pPosCorr_
    (
        IOobject
        (
            lookupOrAddDefault<word>("pPosCorrName", "pPosCorr"),
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
            lookupOrAddDefault<word>("UPosCorrName", "UPosCorr"),
            db.time().timeName(),
            db,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    )
{
    if (debug)
    {
        Info<< dictionary::name() << " = " << *this << endl;
    }
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
    const volScalarField& pnd = cloud.pndcPdfInst();

    volScalarField rhs = rho - pnd;
    // correction factor to make the RHS of the Poisson equation integrate
    // to zero
    dimensionedScalar beta = fvc::domainIntegrate(rhs) / gSum(mesh().V());
    beta.dimensions() /= dimVolume;

    // solve Poisson equation for pPosCorr
    fvScalarMatrix pPosCorrEqn
    (
        fvm::laplacian(pPosCorr_) == C_/sqr(deltaT)*(rhs - beta)
    );
    pPosCorrEqn.setReference(pRefCell_, pRefValue_);
    pPosCorrEqn.relax();
    pPosCorrEqn.solve();

    // time-integrator for correction velocity
    solve(fvm::ddt(rho, UPosCorr_) == -fvc::grad(pPosCorr_));

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
