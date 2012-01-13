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

#include "mcMuradogluPositionCorrection.H"

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcMuradogluPositionCorrection, 0);
    addNamedToRunTimeSelectionTable
    (
        mcPositionCorrection,
        mcMuradogluPositionCorrection,
        mcPositionCorrection,
        Muradoglu
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcMuradogluPositionCorrection::mcMuradogluPositionCorrection
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& dict
)
:
    mcEllipticRelaxationPositionCorrection(db, dict),

    phi_
    (
        IOobject
        (
            lookupOrDefault<word>("phiName", "phiPosCorr"),
            db.time().timeName(),
            db,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::mcMuradogluPositionCorrection::correct(Foam::mcParticleCloud& cloud)
{
    // maximum mean velocity
    dimensionedScalar U0 = max(mag(cloud.Ufv()));
    reduce(U0.value(), maxOp<scalar>());
    // frequency helper variable
    volScalarField omega = c_*U0/L();
    // normalized density difference
    volScalarField QInst =
    (
        cloud.pndcPdfInst() - cloud.rhocPdfInst()
    )/cloud.rhocPdf();

    // solve for mean normalized density difference
    solve
    (
        fvm::ddt(Q_)
        ==
      - fvm::Sp(omega, Q_) + omega*QInst
      + f_*U0*L()*fvm::laplacian(Q_)
    );

    // time-integrator for correction potential
    solve(fvm::ddt(phi_) == b_*sqr(U0)*QInst);

    // indicator variable
    volScalarField zeta = pos(cloud.pndcPdfInst()/cloud.rhocPdfInst() - eps_);

    // note: gradInterpolationConstantTet gives really bad results.
    volVectorField gradPhi = fvc::grad(phi_);
    volVectorField gradQInst = fvc::grad(QInst);
    volVectorField gradQ = fvc::grad(Q_);
    interpolationCellPointFace<vector> gradPhiInterp(gradPhi);
    interpolationCellPointFace<vector> gradQInstInterp(gradQInst);
    interpolationCellPointFace<vector> gradQInterp(gradQ);
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
            gradPhiInterp.interpolate(p, c, f)
          + (a_*U0*l).value()
          * (z != 0 ? gradQInterp : gradQInstInterp).interpolate(p, c, f);
    }
}

// ************************************************************************* //
