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

#include "mcSLMFullVelocityModel.H"

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcSLMFullVelocityModel, 0);
    addNamedToRunTimeSelectionTable
    (
        mcVelocityModel,
        mcSLMFullVelocityModel,
        mcVelocityModel,
        SLMFull
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcSLMFullVelocityModel::mcSLMFullVelocityModel
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& parentDict,
    const Foam::dictionary& dict
)
:
    mcVelocityModel(db, parentDict, dict),
    mesh_(refCast<const fvMesh>(db)),
    pfv_
    (
        db.lookupObject<volScalarField>
        (
            lookupOrDefault<word>("pName", "p", true)
        )
    ),
    p_
    (
        IOobject
        (
            "SLMFullVelocityModel::p",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    diffU_
    (
        IOobject
        (
            "SLMFullVelocityModel::diffU",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimVelocity/dimTime,
        zeroGradientFvPatchScalarField::typeName
    ),
    tetDecomp_(mesh_),
    C0_(lookupOrAddDefault<scalar>("C0", 2.1, true)),
    C1_(lookupOrAddDefault<scalar>("C1", 1.0, true)),
    URelaxTime_
    (
        "URelaxTime",
        dimTime,
        lookupOrAddDefault<scalar>
                ("URelaxTime", db.time().deltaT().value()*10.0, true)
    ),
    kRelaxTime_
    (
        "kRelaxTime",
        dimTime,
        lookupOrAddDefault<scalar>
                ("kRelaxTime", db.time().deltaT().value()*10.0, true)
    )
{
    scalar dt = db.time().deltaT().value();
    // Sanitize input

    C0_ = max(0.0, C0_);
    set("C0", C0_);

    C1_ = max(0.0, min(1.0, C1_));
    set("C1", C1_);

    URelaxTime_.value()   = max(2.0*dt,  URelaxTime_.value());
    set("URelaxTime", URelaxTime_.value());

    kRelaxTime_.value()   = max(2.0*dt,  kRelaxTime_.value());
    set("kRelaxTime", kRelaxTime_.value());

    if (debug)
    {
        Info<< "Sanitized " << dictionary::name() << ":\n"
            << *this << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcSLMFullVelocityModel::setupInternals(const mcParticleCloud& cloud)
{
    const volVectorField& Updf =
        mesh_.lookupObject<volVectorField>("UCloudPDF");

    p_ = pfv_ - 2./3.*cloud.rhocPdf()*cloud.kfv();

    diffU_.internalField() = (cloud.Ufv() - Updf)/URelaxTime_;
    diffU_.correctBoundaryConditions();

    diffk_ = (sqrt(cloud.kfv()/cloud.kcPdf()) - 1.0)/kRelaxTime_;

    rhoInterp_.reset(new interpolationCellPointFace<scalar>(cloud.rhocPdf()));
    UInterp_.reset(new interpolationCellPointFace<vector>(cloud.Ufv()));
    gradPInterp_.reset
    (
        new gradInterpolationConstantTet<scalar>
        (
            tetDecomp_,
            p_
        )
    );
    kInterp_.reset(new interpolationCellPointFace<scalar>(cloud.kfv()));
    diffUInterp_.reset(new interpolationCellPointFace<vector>(diffU_));
    kcPdfInterp_.reset(new interpolationCellPointFace<scalar>(cloud.kcPdf()));
}


void Foam::mcSLMFullVelocityModel::correct
(
    Foam::mcParticleCloud& cloud,
    Foam::mcParticle& p,
    bool prepare
)
{
    if (prepare)
    {
        setupInternals(cloud);
    }

    scalar deltaT = p.eta()*mesh_.time().deltaT().value();

    const point& pos = p.position();
    label c = p.cell();
    label f = p.face();

    // fluid quantities @ particle position
    vector UFap = UInterp_().interpolate(pos, c, f);
    vector gradPFap = gradPInterp_().interpolate(pos, c, f);
    scalar kFap = kInterp_().interpolate(pos, c, f);
    vector diffUap = diffUInterp_().interpolate(pos, c, f);

    // Note: it would be the best if UInterp was interpolating velocities based
    // on face fluxes instead of cell center values. Will implement later.

    const vector xi =
        vector
        (
            cloud.random().GaussNormal(),
            cloud.random().GaussNormal(),
            cloud.random().GaussNormal()
        );

    const scalar A = -(0.5*C1_ + 0.75*C0_)*p.Omega();
    const vector B = -(gradPFap/p.rho() + A*UFap);
    const vector deltaU =
        (B + A*p.UParticle())*deltaT + sqrt(C0_*kFap*p.Omega()*deltaT)*xi;

    p.UParticle() += deltaU + 0.5*A*deltaU*deltaT
        // Correct mean velocity
        + diffUap*deltaT
        // Scale to ensure consistency on TKE (using interpolated
        // kFap/kPdfap gives very bad results)
        + (p.UParticle() - UFap)*diffk_[c]*deltaT;
}

// ************************************************************************* //
