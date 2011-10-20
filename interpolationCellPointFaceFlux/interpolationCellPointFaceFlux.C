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

Description

\*---------------------------------------------------------------------------*/

#include "interpolationCellPointFaceFlux.H"
#include "volFields.H"
#include "polyMesh.H"
#include "volPointInterpolation.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  // Private member functions

  tmp<surfaceVectorField> interpolationCellPointFaceFlux::makeFaceFeild
  (
   const volVectorField& psi,
   const word& fluxName
   )
  {
    tmp<surfaceVectorField> tresult = linearInterpolate(psi);
    surfaceVectorField& result = tresult();

    const surfaceScalarField& phi = psi.mesh().lookupObject<surfaceScalarField>(fluxName);

    surfaceVectorField n = mesh.Sf()/mesh.magSf();

    result -= n*(n & result);
    result += phi/mesh.magSf()*n;

    return tresult;
  }

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

interpolationCellPointFaceFlux::interpolationCellPointFaceFlux
(
    const volVectorField& psi
)
:
    interpolationCellPointFace<vector>
    (
        psi,
        volPointInterpolation::New(psi.mesh()).interpolate(psi),
        makeFaceField(psi, "phi")
    )
{
  Info << " calling new constructor. " << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
