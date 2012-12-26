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

#include "mcEllipticRelaxationPositionCorrection.H"

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "interpolation.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcEllipticRelaxationPositionCorrection, 0);
    addNamedToRunTimeSelectionTable
    (
        mcPositionCorrection,
        mcEllipticRelaxationPositionCorrection,
        mcPositionCorrection,
        ellipticRelaxation
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcEllipticRelaxationPositionCorrection::
mcEllipticRelaxationPositionCorrection
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcPositionCorrection(cloud, db, subDictName),

    Q_
    (
        IOobject
        (
            thermoDict().lookupOrDefault<word>("QName", "QPosCorr"),
            db.time().timeName(),
            db,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        cloud.mesh()
    ),

    QInst_
    (
        IOobject
        (
            "mcEllipticRelaxationPositionCorrection::QInst",
            db.time().timeName(),
            db,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud.mesh(),
        dimless,
        zeroGradientFvPatchScalarField::typeName
    ),

    gradQInst_
    (
        IOobject
        (
            "mcEllipticRelaxationPositionCorrection::grad(QInst)",
            db.time().timeName(),
            db,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud.mesh(),
        dimless/dimLength,
        zeroGradientFvPatchScalarField::typeName
    ),

    gradQ_
    (
        IOobject
        (
            "mcEllipticRelaxationPositionCorrection::grad(Q)",
            db.time().timeName(),
            db,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud.mesh(),
        dimless/dimLength,
        zeroGradientFvPatchScalarField::typeName
    ),

    zeta_
    (
        IOobject
        (
            "mcEllipticRelaxationPositionCorrection::zeta",
            db.time().timeName(),
            db,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud.mesh(),
        dimless,
        zeroGradientFvPatchScalarField::typeName
    ),

    U0_("U0", dimVelocity, 0.),
    eps_("eps", dimless, 0.),
    c_("c", dimless, 0.),
    f_("f", dimless, 0.),
    b_("b", dimless, 0.),
    a_("a", dimless, 0.),
    LInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solutionDict().interpolationScheme(L().name()),
            L()
        )
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcEllipticRelaxationPositionCorrection::readCoeffs()
{
#if FOAM_HEX_VERSION < 0x200
    using mathematicalConstant::pi;
#else
    using constant::mathematical::pi;
#endif
    dimensionedScalar kf
    (
        "kf",
        dimless,
        solutionDict().lookupOrDefault<scalar>("kf", 3.)
    );
    dimensionedScalar kb
    (
        "kb",
        dimDensity,
        solutionDict().lookupOrDefault<scalar>("kb", 8.)
    );
    dimensionedScalar N
    (
        "N",
        dimless,
        solutionDict().lookupOrDefault<scalar>("N", 20.)
    );
    dimensionedScalar CFL
    (
        "CFL",
        dimless,
        solutionDict().lookupOrDefault<scalar>("CFL", 0.4)
    );
    eps_ = solutionDict().lookupOrDefault<scalar>("eps", 0.25);
    c_ = solutionDict().lookupOrDefault<scalar>
        (
            "c",
            1./(pi*CFL*N).value()
        );
    f_ = solutionDict().lookupOrDefault<scalar>("f", (kf*c_).value());
    b_ = solutionDict().lookupOrDefault<scalar>("b", (kb*sqr(f_)).value());
    a_ = solutionDict().lookupOrDefault<scalar>
        (
            "a",
            f_.value()*(1. + (b_/sqr(c_)).value())
        );
}

void Foam::mcEllipticRelaxationPositionCorrection::updateInternals()
{
    readCoeffs();
    U0_ = max(mag(cloud().Ufv()));

    QInst_ = (cloud().pndcPdfInst() - cloud().rhocPdfInst())/cloud().rhocPdf();

    // elliptic relaxation
    fvScalarMatrix QEqn
    (
        fvm::Sp(1, Q_)
      - f_*sqr(L())/c_*fvm::laplacian(Q_)
        ==
        QInst_
    );

    QEqn.setReference(0, 0.0);
    QEqn.solve();

    zeta_ = pos(cloud().pndcPdfInst()/cloud().rhocPdfInst() - eps_);
    gradQ_ = fvc::grad(Q_);
    gradQInst_ = fvc::grad(QInst_);
    gradQInterp_ = interpolation<vector>::New
    (
        cloud().solutionDict().subDict("interpolationSchemes"),
        gradQ_
    );
    gradQInstInterp_  = interpolation<vector>::New
    (
        cloud().solutionDict().subDict("interpolationSchemes"),
        gradQInst_
    );
    zetaInterp_  = interpolation<scalar>::New
    (
        cloud().solutionDict().subDict("interpolationSchemes"),
        zeta_
    );
}


void Foam::mcEllipticRelaxationPositionCorrection::correct(mcParticle& part)
{
    const point& p = part.position();
    label c = part.cell();
    label f = part.face();
    scalar z = zetaInterp_().interpolate(p, c, f);
    dimensionedScalar l = LInterp_().interpolate(p, c, f);
    part.Ucorrection() -=
      + (a_*U0_*l).value()
      * (z != 0 ? gradQInterp_ : gradQInstInterp_)().interpolate(p, c, f);
}

// ************************************************************************* //
