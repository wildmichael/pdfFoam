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

#include "mcPositionCorrection.H"
#include "mcParticleCloud.H"

#include "fvCFD.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcPositionCorrection, 0);
    defineRunTimeSelectionTable(mcPositionCorrection, mcPositionCorrection);

} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcPositionCorrection::mcPositionCorrection
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcModel(cloud, db, subDictName)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mcPositionCorrection> Foam::mcPositionCorrection::New
(
    mcParticleCloud& cloud,
    const objectRegistry& db
)
{
    word posCorrType(cloud.thermoDict().lookup("positionCorrection"));
    word sd = posCorrType + "PositionCorrectionCoeffs";

    // If set to "false", "no", "n", "off" or "none", disable
#if FOAM_HEX_VERSION < 0x200
    Switch::switchType enablePosCorr = Switch::asEnum(posCorrType, true);
    if (enablePosCorr != Switch::INVALID && !Switch::asBool(enablePosCorr))
#else
    Switch enablePosCorr = Switch(posCorrType, true);
    if (enablePosCorr.valid() && !bool(enablePosCorr))
#endif
    {
        return autoPtr<mcPositionCorrection>
        (
            new mcPositionCorrection(cloud, db, sd)
        );
    }

    mcPositionCorrectionConstructorTable::iterator cstrIter =
        mcPositionCorrectionConstructorTablePtr_->find(posCorrType);

    if (cstrIter == mcPositionCorrectionConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mcPositionCorrection::New(const fvMesh&, const dictionary&)"
        )   << "Unknown mcPositionCorrection type " << posCorrType << nl << nl
            << "Valid mcPositionCorrection types are :" << nl
            << mcPositionCorrectionConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<mcPositionCorrection>
    (
        cstrIter()(cloud, db, sd)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcPositionCorrection::makeL()
{
#if FOAM_HEX_VERSION < 0x200
    using mathematicalConstant::pi;
#else
    using constant::mathematical::pi;
#endif
    const fvMesh& mesh = cloud().mesh();
    if (L_.valid())
    {
        return;
    }

    L_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "mcPositionCorrecton::L",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("0", dimLength, 0.),
            calculatedFvPatchScalarField::typeName
        )
    );

    const surfaceScalarField& invdx =
        mesh.surfaceInterpolation::deltaCoeffs();
    forAll(mesh.cells(), cellI)
    {
        const cell& c = mesh.cells()[cellI];
        scalar sumDx = 0.;
        label n = 0;
        forAll(c, cellFaceI)
        {
            label faceI = c[cellFaceI];
            scalar coeff;
            if (mesh.isInternalFace(faceI))
            {
                coeff = invdx[faceI];
            }
            else
            {
                label patchI = mesh.boundaryMesh().whichPatch(faceI);
                if (!invdx.boundaryField()[patchI].size())
                {
                    // skip empty boundaries
                    continue;
                }
                label start = mesh.boundaryMesh()[patchI].start();
                coeff = invdx.boundaryField()[patchI][faceI-start]/2.;
            }
            sumDx += 1./coeff;
            ++n;
        }
        L_()[cellI] = sumDx/(n*pi);
    }
    forAll(L_().boundaryField(), patchI)
    {
        L_().boundaryField()[patchI] =
            2./(invdx.boundaryField()[patchI]*pi);
    }
}


void Foam::mcPositionCorrection::correct(mcParticle& p)
{}

// ************************************************************************* //
