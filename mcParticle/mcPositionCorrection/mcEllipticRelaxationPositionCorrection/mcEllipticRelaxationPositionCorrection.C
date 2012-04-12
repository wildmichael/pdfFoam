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

#include "mcEllipticRelaxationPositionCorrection.H"

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
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
    const Foam::objectRegistry& db,
    const Foam::dictionary& parentDict,
    const Foam::dictionary& dict
)
:
    mcPositionCorrection(db, parentDict, dict),

    Q_
    (
        IOobject
        (
            lookupOrDefault<word>("QName", "QPosCorr"),
            db.time().timeName(),
            db,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),

    kf_
    (
        "kf",
        dimless,
        lookupOrAddDefault<scalar>("kf", 3.)
    ),
    kb_
    (
        "kb",
        dimDensity,
        lookupOrAddDefault<scalar>("kb", 8.)
    ),
    N_
    (
        "N",
        dimless,
        lookupOrAddDefault<scalar>("N", 20.)
    ),
    CFL_
    (
        "CFL",
        dimless,
        lookupOrAddDefault<scalar>("CFL", 0.4)
    ),
    eps_
    (
        "eps",
        dimless,
        lookupOrAddDefault<scalar>("eps", 0.25)
    ),
    c_
    (
        "c",
        dimless,
        lookupOrAddDefault<scalar>
        (
            "c",
            1./(mathematicalConstant::pi*CFL_*N_).value()
        )
    ),
    f_
    (
        "f",
        dimless,
        lookupOrAddDefault<scalar>
        (
            "f",
            (kf_*c_).value()
        )
    ),
    b_
    (
        "b",
        dimless,
        lookupOrAddDefault<scalar>
        (
            "b",
            (kb_*sqr(f_)).value()
        )
    ),
    a_
    (
        "a",
        dimless,
        lookupOrAddDefault<scalar>
        (
            "a",
            f_.value()*(1. + (b_/sqr(c_)).value())
        )
    ),
    LInterp_(L())
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::mcEllipticRelaxationPositionCorrection::correct
(
    Foam::mcParticleCloud& cloud
)
{
    // maximum mean velocity
    dimensionedScalar U0 = max(mag(cloud.Ufv()));

    // normalized density difference
    volScalarField QInst =
    (
        cloud.pndcPdfInst() - cloud.rhocPdfInst()
    )/cloud.rhocPdf();

    // elliptic relaxation
    fvScalarMatrix QEqn
    (
        fvm::Sp(1, Q_)
      - f_*sqr(L())/c_*fvm::laplacian(Q_)
        ==
        QInst
    );

    QEqn.setReference(0, 0.0);
    QEqn.solve();

    // indicator variable
    volScalarField zeta = pos(cloud.pndcPdfInst()/cloud.rhocPdfInst() - eps_);

    // TODO try gradInterpolationConstantTet
    volVectorField gradQ = fvc::grad(Q_);
    volVectorField gradQInst = fvc::grad(QInst);
    interpolationCellPointFace<vector> gradQInterp(gradQ);
    interpolationCellPointFace<vector> gradQInstInterp(gradQInst);
    interpolationCellPointFace<scalar> zetaInterp(zeta);

    // apply
    forAllIter(mcParticleCloud, cloud, pIter)
    {
        mcParticle& part = pIter();
        const point& p = part.position();
        label c = part.cell();
        label f = part.face();
        scalar z = zetaInterp.interpolate(p, c, f);
        dimensionedScalar l = LInterp_.interpolate(p, c, f);
        part.Ucorrection() -=
          + (a_*U0*l).value()
          * (z != 0 ? gradQInterp : gradQInstInterp).interpolate(p, c, f);
    }
}

// ************************************************************************* //
