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

#include "mcLimitedSimplePositionCorrection.H"

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcLimitedSimplePositionCorrection, 0);
    addNamedToRunTimeSelectionTable
    (
        mcPositionCorrection,
        mcLimitedSimplePositionCorrection,
        mcPositionCorrection,
        limitedSimple
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcLimitedSimplePositionCorrection::
mcLimitedSimplePositionCorrection
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& parentDict,
    const Foam::dictionary& dict
)
:
    mcPositionCorrection(db, parentDict, dict),

    C_
    (
        "C",
        dimless,
        lookupOrAddDefault<scalar>("C", 1.)
    ),

    CFL_
    (
        "CFL",
        dimless,
        lookupOrAddDefault<scalar>("CFL", 0.5)
    ),

    Cd_
    (
        "Cd",
        dimless,
        lookupOrAddDefault<scalar>("Cd", 0.)
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::mcLimitedSimplePositionCorrection::correct
(
    Foam::mcParticleCloud& cloud
)
{
    scalar dt = mesh().time().deltaT().value();
    volScalarField phi
    (
        IOobject
        (
            "phiPosCorr",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimless, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    );
    phi.internalField() = (cloud.pndcPdf() - cloud.rhocPdf())/cloud.rhocPdf();
    phi.correctBoundaryConditions();

    // limiter
    scalar phiLim = 0.;
    forAll(phi, i)
    {
        phiLim = max(phiLim, fabs(phi[i]));
    }
    reduce(phiLim, maxOp<scalar>());
    phiLim = tanh(10.*phiLim);

    // CFL
    volVectorField UPosCorr = -fvc::grad(phi);
    scalar Co = 0;
    forAll(mesh().cells(), cellI)
    {
        const cell& c = mesh().cells()[cellI];
        forAll(c, cellFaceI)
        {
            label faceI = c[cellFaceI];
            vector coeff;
            if (mesh().isInternalFace(faceI))
            {
                coeff = cloud.CourantCoeffs()[faceI];
            }
            else
            {
                label patchI = mesh().boundaryMesh().whichPatch(faceI);
                if (!cloud.CourantCoeffs().boundaryField()[patchI].size())
                {
                    // skip empty boundaries
                    continue;
                }
                label start = mesh().boundaryMesh()[patchI].start();
                coeff = cloud.CourantCoeffs().boundaryField()[patchI][faceI-start];
            }
            Co = max(Co, fabs(dt*coeff&UPosCorr[cellI]));
        }
    }
    reduce(Co, maxOp<scalar>());
    dimensionedScalar corr = 0.;
    if (Co > SMALL)
    {
        corr = CFL_/Co*C_*phiLim;
    }

    UPosCorr *= corr;
    interpolationCellPointFace<vector> UPosCorrInterp(UPosCorr);

    // TODO try gradInterpolationConstantTet
    //phi *= corr;
    //gradinterpolationConstantTet<scalar> UPosCorrInterp(phi);

    // apply
    forAllIter(mcParticleCloud, cloud, pIter)
    {
        mcParticle& part = pIter();
        const point& p = part.position();
        label c = part.cell();
        label f = part.face();
        vector Uc = UPosCorrInterp.interpolate(p, c, f);
        vector xi
            (
                cloud.random().GaussNormal(),
                cloud.random().GaussNormal(),
                cloud.random().GaussNormal()
            );
        Uc += (Cd_*sqrt(mag(Uc)*dt)*xi/dt).value();
        part.Ucorrection() += Uc;
    }
}

// ************************************************************************* //
