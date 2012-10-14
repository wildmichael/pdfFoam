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
#include "interpolation.H"
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

Foam::mcLimitedSimplePositionCorrection::mcLimitedSimplePositionCorrection
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcPositionCorrection(cloud, db, subDictName),

    UPosCorr_
    (
        IOobject
        (
            thermoDict().lookupOrDefault<word>("UPosCorrName", "UPosCorr"),
            db.time().timeName(),
            db,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cloud.mesh(),
        dimVelocity,
        cloud.rhocPdf().boundaryField().types()
    )
{
    // Set fixedValue boundaries of UPosCorr to vector::zero
    forAll(UPosCorr_.boundaryField(), patchi)
    {
        fvPatchField<vector>& fpf = UPosCorr_.boundaryField()[patchi];
        if (isA<fixedValueFvPatchVectorField>(fpf))
        {
            fpf.vectorField::operator=(vector::zero);
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcLimitedSimplePositionCorrection::updateInternals()
{
    dimensionedScalar C
    (
        "C",
        dimless,
        solutionDict().lookupOrDefault<scalar>("C", 1e-2)
    );

    dimensionedScalar CFL
    (
        "CFL",
        dimless,
        solutionDict().lookupOrDefault<scalar>("CFL", 0.5)
    );

    static const dimensionedScalar one("one", dimLength*dimVelocity, 1.);

    const fvMesh& mesh = cloud().mesh();
    scalar dt = cloud().deltaT().value();
    volScalarField phi
    (
        IOobject
        (
            "phiPosCorr",
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        0.*one,
        zeroGradientFvPatchField<scalar>::typeName
    );
    phi.internalField() =
        one*(cloud().pndcPdf() - cloud().rhocPdf())/cloud().rhocPdf();
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
    UPosCorr_ = -fvc::grad(phi);
    scalar Co = 0;
    forAll(mesh.cells(), cellI)
    {
        const cell& c = mesh.cells()[cellI];
        forAll(c, cellFaceI)
        {
            label faceI = c[cellFaceI];
            vector coeff;
            if (mesh.isInternalFace(faceI))
            {
                coeff = cloud().CourantCoeffs()[faceI];
            }
            else
            {
                label patchI = mesh.boundaryMesh().whichPatch(faceI);
                if (!cloud().CourantCoeffs().boundaryField()[patchI].size())
                {
                    // skip empty boundaries
                    continue;
                }
                label start = mesh.boundaryMesh()[patchI].start();
                coeff =
                    (
                        cloud().CourantCoeffs()
                            .boundaryField()[patchI][faceI-start]
                    );
            }
            Co = max(Co, fabs(dt*coeff&UPosCorr_[cellI]));
        }
    }
    reduce(Co, maxOp<scalar>());
    dimensionedScalar corr = 0.;
    if (Co > SMALL)
    {
        corr = CFL/Co*C*phiLim;
    }

    UPosCorr_ *= corr;
    UPosCorr_.correctBoundaryConditions();
    UPosCorrInterp_ = interpolation<vector>::New
    (
        cloud().solutionDict().interpolationScheme(UPosCorr_.name()),
        UPosCorr_
    );

    // TODO try gradInterpolationConstantTet
    //phi *= corr;
    //gradinterpolationConstantTet<scalar> UPosCorrInterp(phi);
}


void Foam::mcLimitedSimplePositionCorrection::correct(mcParticle& part)
{
    const point& p = part.position();
    label c = part.cell();
    label f = part.face();
    part.Ucorrection() += UPosCorrInterp_().interpolate(p, c, f);
    if (debug)
    {
        scalar UPosCorrMax = gMax(mag(UPosCorr_)());
        scalar UCorrMag = mag(part.Ucorrection());
        if (UCorrMag > 1.5*UPosCorrMax)
        {
            WarningIn
            (
                "mcLimitedSimplePositionCorrection::correct(mcParticle&)"
            )   << "Interpolation of UPosCorr to location " << p
                << " yielded a correction velocity magnitude of "
                << UCorrMag << ", where the maximum magnitude in the field is "
                << UPosCorrMax << endl;
        }
    }
}

// ************************************************************************* //
