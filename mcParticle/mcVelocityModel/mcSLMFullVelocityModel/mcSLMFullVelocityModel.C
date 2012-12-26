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

\*---------------------------------------------------------------------------*/

#include "mcSLMFullVelocityModel.H"

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "interpolation.H"
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
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcVelocityModel(cloud, db, subDictName),
    pfv_
    (
        db.lookupObject<volScalarField>
        (
            thermoDict().lookupOrDefault<word>("pName", "p")
        )
    ),
    p_
    (
        IOobject
        (
            "SLMFullVelocityModel::p",
            cloud.mesh().time().timeName(),
            cloud.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud.mesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    diffU_
    (
        IOobject
        (
            "SLMFullVelocityModel::diffU",
            cloud.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud.mesh(),
        dimVelocity/dimTime,
        zeroGradientFvPatchScalarField::typeName
    ),
    tetDecomp_(cloud.mesh()),
    C0_(thermoDict().lookupOrDefault<scalar>("C0", 2.1)),
    C1_(thermoDict().lookupOrDefault<scalar>("C1", 1.0))
{
    // Sanitize input
    C0_ = max(0.0, C0_);
    C1_ = max(0.0, min(1.0, C1_));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcSLMFullVelocityModel::updateInternals()
{
    const volVectorField& Updf =
        cloud().mesh().lookupObject<volVectorField>("UCloudPDF");
    const mcSolution& solDict = cloud().solutionDict();

    p_ = pfv_ - 2./3.*cloud().rhocPdf()*cloud().kfv();

    diffU_.internalField() = (cloud().Ufv() - Updf)/solDict.relaxationTime("U");
    diffU_.correctBoundaryConditions();

    diffk_ = (sqrt(cloud().kfv()/cloud().kcPdf()) - 1.0)
            /solDict.relaxationTime("k");

    const volVectorField& Ufv = cloud().Ufv();
    UInterp_ = interpolation<vector>::New
    (
        solDict.interpolationScheme(Ufv.name()),
        Ufv
    );
    gradPInterp_.reset
    (
        new gradInterpolationConstantTet<scalar>
        (
            tetDecomp_,
            p_
        )
    );
    const volScalarField& kfv = cloud().kfv();
    kInterp_ = interpolation<scalar>::New
    (
        solDict.interpolationScheme(kfv.name()),
        kfv
    );
    diffUInterp_ = interpolation<vector>::New
    (
        solDict.interpolationScheme(diffU_.name()),
        diffU_
    );
}


void Foam::mcSLMFullVelocityModel::correct(mcParticle& p)
{
    const scalar& deltaTg = cloud().deltaT().value();
    scalar deltaT = p.eta()*deltaTg;

    const point& pos = p.position();
    label c = p.cell();
    label f = p.face();

    // fluid quantities @ particle position
    vector UFap = UInterp_().interpolate(pos, c, f);
    vector gradPFap = gradPInterp_().interpolate(pos, c, f);
    const scalar& kMin = cloud().solutionDict().kMin().value();
    scalar kFap = max(kInterp_().interpolate(pos, c, f), kMin);
    vector diffUap = diffUInterp_().interpolate(pos, c, f);

    // Note: it would be the best if UInterp was interpolating velocities based
    // on face fluxes instead of cell center values. Will implement later.

    const vector xi =
        vector
        (
            cloud().random().GaussNormal(),
            cloud().random().GaussNormal(),
            cloud().random().GaussNormal()
        );

    const scalar A = -(0.5*C1_ + 0.75*C0_)*p.Omega();
    const vector B = -(gradPFap/p.rho() + A*UFap);
    const vector deltaU =
        (B + A*p.UParticle())*deltaT + sqrt(C0_*kFap*p.Omega()*deltaT)*xi;

    p.UParticle() += deltaU + 0.5*A*deltaU*deltaT
        // Correct mean velocity
        + diffUap*deltaTg
        // Scale to ensure consistency on TKE (using interpolated
        // kFap/kPdfap gives very bad results)
        + (p.UParticle() - UFap)*diffk_[c]*deltaTg;
    p.Co() = max(p.Co(), p.Omega());
}

// ************************************************************************* //
