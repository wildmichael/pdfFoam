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

#include "mcPositionCorrection.H"

#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcPositionCorrection, 0);
    defineRunTimeSelectionTable(mcPositionCorrection, mcPositionCorrection);

} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcPositionCorrection::mcPositionCorrection
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& parentDict,
    const Foam::dictionary& mcPositionCorrectionDict
)
:
    mcModel(db, parentDict, mcPositionCorrectionDict),
    mesh_(refCast<const fvMesh>(db))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mcPositionCorrection> Foam::mcPositionCorrection::New
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& dict
)
{
    word posCorrType(dict.lookup("positionCorrection"));

    // If set to "false", "no", "n", "off" or "none", disable
    Switch::switchType enablePosCorr = Switch::asEnum(posCorrType, true);
    if (enablePosCorr != Switch::INVALID && !Switch::asBool(enablePosCorr))
    {
        return autoPtr<mcPositionCorrection>
        (
            new mcPositionCorrection(db, dict, dictionary::null)
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
        cstrIter()
        (
            db,
            dict,
            dict.subOrEmptyDict(posCorrType+"PositionCorrectionCoeffs")
        )
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcPositionCorrection::makeL()
{

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
                "LPosCorr_",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimLength, 0.),
            calculatedFvPatchScalarField::typeName
        )
    );

    const surfaceScalarField& invdx =
        mesh_.surfaceInterpolation::deltaCoeffs();
    forAll(mesh_.cells(), cellI)
    {
        const cell& c = mesh_.cells()[cellI];
        scalar sumDx = 0.;
        label n = 0;
        forAll(c, cellFaceI)
        {
            label faceI = c[cellFaceI];
            scalar coeff;
            if (mesh_.isInternalFace(faceI))
            {
                coeff = invdx[faceI];
            }
            else
            {
                label patchI = mesh_.boundaryMesh().whichPatch(faceI);
                if (!invdx.boundaryField()[patchI].size())
                {
                    // skip empty boundaries
                    continue;
                }
                label start = mesh_.boundaryMesh()[patchI].start();
                coeff = invdx.boundaryField()[patchI][faceI-start]/2.;
            }
            sumDx += 1./coeff;
            ++n;
        }
        L_()[cellI] = sumDx/(n*mathematicalConstant::pi);
    }
    forAll(L_().boundaryField(), patchI)
    {
        L_().boundaryField()[patchI] =
            2./(invdx.boundaryField()[patchI]*mathematicalConstant::pi);
    }
}


void Foam::mcPositionCorrection::correct
(
    Foam::mcParticleCloud& cloud
)
{}

// ************************************************************************* //
