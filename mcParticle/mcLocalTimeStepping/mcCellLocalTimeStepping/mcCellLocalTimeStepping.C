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

#include "mcCellLocalTimeStepping.H"

#include "addToRunTimeSelectionTable.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "fvCFD.H"
#include "interpolation.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcCellLocalTimeStepping, 0);
    addNamedToRunTimeSelectionTable
    (
        mcLocalTimeStepping,
        mcCellLocalTimeStepping,
        mcLocalTimeStepping,
        cell
    );

} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcCellLocalTimeStepping::mcCellLocalTimeStepping
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcLocalTimeStepping(cloud, db, subDictName),
    eta_
    (
        IOobject
        (
            "mcCellLocaltimeStepping::eta",
            db.time().timeName(),
            db,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        refCast<const fvMesh>(db),
        dimless,
        zeroGradientFvPatchScalarField::typeName
    )
{
    updateInternals();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcCellLocalTimeStepping::updateInternals()
{
    // Maximum local to global time step ratio
    scalar etaMax =
        max(solutionDict().lookupOrDefault<scalar>("upperBound", 10), 1);
    const fvMesh& mesh = eta_.mesh();
    scalar deltaT = cloud().deltaT().value();
    const volVectorField& U = cloud().Ufv();
    const compressible::turbulenceModel& tm = cloud().turbulenceModel();
    tmp<volSymmTensorField> tR = tm.R();
    const symmTensorField& RInt = tR().internalField();
    const surfaceVectorField& CourantCoeffs = cloud().CourantCoeffs();
    scalar k0 = SMALL;
    if (isA<compressible::RASModel>(tm))
    {
        k0 = refCast<const compressible::RASModel>(tm).k0().value();
    }

    // Use clipping in case there are negative entries in R due to divergence
    // errors
    vectorField urms2 =
    (
        vector(2, 0, 0) * sqrt(max(k0, RInt.component(symmTensor::XX)))
      + vector(0, 2, 0) * sqrt(max(k0, RInt.component(symmTensor::YY)))
      + vector(0, 0, 2) * sqrt(max(k0, RInt.component(symmTensor::ZZ)))
    );

    forAll(eta_, cellI)
    {
        scalar dtc = GREAT;
        const cell& c = mesh.cells()[cellI];

        forAll(c, cellFaceI)
        {
           label faceI = c[cellFaceI];
           vector coeff;
           if (mesh.isInternalFace(faceI))
           {
               coeff = CourantCoeffs[faceI];
           }
           else
           {
               label patchI = mesh.boundaryMesh().whichPatch(faceI);
               if (!CourantCoeffs.boundaryField()[patchI].size())
               {
                   // skip empty boundaries
                   continue;
               }
               label start = mesh.boundaryMesh()[patchI].start();
               coeff = CourantCoeffs.boundaryField()[patchI][faceI-start];
           }
           dtc =
               min
               (
                   dtc,
                   1./
                   (
                       fabs(coeff&U[cellI])
                     + fabs(coeff&urms2[cellI])
                     + SMALL
                   )
               );
        }
        eta_[cellI] = dtc/deltaT;
    }
    // scale eta field
    scalarField& etaInt = eta_.internalField();
    scalar etaMin = gMin(etaInt);
    etaInt -= etaMin;
    etaInt = (etaMax-1)/gMax(etaInt)*etaInt + 1;
    eta_.correctBoundaryConditions();
    etaInterp_ = interpolation<scalar>::New
    (
        cloud().solutionDict().interpolationScheme(eta_.name()),
        eta_
    );
}


void Foam::mcCellLocalTimeStepping::correct(mcParticle& p)
{
    if (!etaInterp_.valid())
    {
        FatalErrorIn("mcCellLocalTimeStepping::correct"
            "(mcParticleCloud&,mcParticle&,bool)")
            << "Interpolator not initialized"
            << exit(FatalError);
    }
    p.eta() = etaInterp_().interpolate(p.position(), p.cell(), p.face());
}

// ************************************************************************* //
