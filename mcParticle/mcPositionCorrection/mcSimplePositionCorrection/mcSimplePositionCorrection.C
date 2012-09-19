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

#include "mcSimplePositionCorrection.H"

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "interpolation.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcSimplePositionCorrection, 0);
    addNamedToRunTimeSelectionTable
    (
        mcPositionCorrection,
        mcSimplePositionCorrection,
        mcPositionCorrection,
        simple
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcSimplePositionCorrection::
mcSimplePositionCorrection
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcPositionCorrection(cloud, db, subDictName),

    Ainv_
    (
        IOobject
        (
            "mcSimplePositionCorrection::A",
            db.time().constant(),
            db,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud.mesh(),
        dimensionedVector("0", dimless/dimArea, vector::zero)
    ),

    phi_
    (
        IOobject
        (
            "mcSimplePositionCorrection::phi_",
            db.time().timeName(),
            db,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud.mesh(),
        dimless,
        zeroGradientFvPatchScalarField::typeName
    ),

    gradPhi_
    (
        IOobject
        (
            "mcSimplePositionCorrection::grad(phi)",
            db.time().timeName(),
            db,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud.mesh(),
        dimless/dimLength,
        zeroGradientFvPatchScalarField::typeName
    )
{
    const fvMesh& m = cloud.mesh();
    const pointField& points = m.points();
    const faceList& f = m.faces();
    const cellList& cf = m.cells();

    forAll(cf, cellI)
    {
        labelList cellVertices = cf[cellI].labels(f);

        vector bbmax = -GREAT*vector::one;
        vector bbmin = GREAT*vector::one;

        forAll (cellVertices, vertexI)
        {
            bbmax = max(bbmax, points[cellVertices[vertexI]]);
            bbmin = min(bbmin, points[cellVertices[vertexI]]);
        }

        scalar dx = bbmax.x()-bbmin.x();
        scalar dy = bbmax.y()-bbmin.y();
        scalar dz = bbmax.z()-bbmin.z();

        Ainv_[cellI] = vector(1./dy*dz, 1./dx*dz, 1./dx*dy);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcSimplePositionCorrection::updateInternals()
{
    dimensionedScalar C
    (
        "C",
        dimless,
        solutionDict().lookupOrDefault<scalar>("C", 1e-2)
    );
    const scalarField& pndInt = cloud().pndcPdf().dimensionedInternalField();
    const scalarField& rhoInt = cloud().rhocPdf().dimensionedInternalField();
    phi_.internalField() = pndInt/rhoInt - 1.;
    phi_.correctBoundaryConditions();

    // TODO try gradInterpolationConstantTet
    gradPhi_ = C*mag(cloud().Ufv())*fvc::grad(phi_);
    gradPhiInterp_ = interpolation<vector>::New
    (
        cloud().solutionDict().interpolationScheme(gradPhi_.name()),
        gradPhi_
    );
}


void Foam::mcSimplePositionCorrection::correct(Foam::mcParticle& part)
{
    const point& p = part.position();
    label c = part.cell();
    label f = part.face();
    part.Ucorrection() -=
        cmptMultiply(Ainv_[c], gradPhiInterp_().interpolate(p, c, f));
}

// ************************************************************************* //
