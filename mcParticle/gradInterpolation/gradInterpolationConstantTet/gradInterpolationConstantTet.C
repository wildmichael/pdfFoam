/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "gradInterpolationConstantTet.H"
#include "volFields.H"
#include "polyMesh.H"
#include "volPointInterpolation.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
gradInterpolationConstantTet<Type>::gradInterpolationConstantTet
(
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
:
    gradInterpolation<Type>(psi),
    psip_(volPointInterpolation::New(psi.mesh()).interpolate(psi)),
    psis_(linearInterpolate(psi)),
    tetDecomp_(new tetFacePointCellDecomposition<richTetPointRef>(psi.mesh()))
{}

template<class Type>
gradInterpolationConstantTet<Type>::gradInterpolationConstantTet
(
    const tetFacePointCellDecomposition<richTetPointRef>& tetDecomp,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
:
    gradInterpolation<Type>(psi),
    psip_(volPointInterpolation::New(psi.mesh()).interpolate(psi)),
    psis_(linearInterpolate(psi)),
    tetDecomp_(tetDecomp)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::outerProduct<Foam::vector, Type>::type
gradInterpolationConstantTet<Type>::interpolate
(
    const vector& position,
    const label cellI,
    const label faceI
) const
{
    typename Foam::outerProduct<Foam::vector, Type>::type grad;
    label tetI = tetDecomp_().find(position, cellI);
    if (tetI > -1)
    {
        label gFaceI = tetDecomp_().tetrahedronFace()[tetI];
        const face& f = this->pMeshFaces_[gFaceI];
        label ptBI = f[tetDecomp_().tetrahedronPoints()[tetI].first()];
        label ptCI = f[tetDecomp_().tetrahedronPoints()[tetI].second()];
        const vectorField& gradNi =
            tetDecomp_().tetrahedra()[tetI].gradNi();
        grad =
        (
            gradNi[3]*this->psi_[cellI]
          + gradNi[1]*psip_[ptBI]
          + gradNi[2]*psip_[ptCI]
        );
        if (this->pMesh_.isInternalFace(gFaceI))
        {
            grad += gradNi[0]*psis_[gFaceI];
        }
        else
        {
            label patchI = this->pMesh_.boundaryMesh().whichPatch(gFaceI);
            if (psis_.boundaryField()[patchI].size())
            {
                label start = this->pMesh_.boundaryMesh()[patchI].start();
                grad += gradNi[0]*psis_.boundaryField()[patchI][gFaceI-start];
            }
            else
            {
                grad += gradNi[0]*this->psi_[cellI];
            }
        }
    }
    else
    {
        label realCellI = this->pMesh_.findCell(position);
        if (realCellI > -1)
        {
            FatalErrorIn("gradInterpolationConstantTet<Type>::interpolate"
                        "("
                            "const vector&, "
                            "const label, "
                            "const label"
                        ") const")
                << "Failed to find tetrahedron containing point " << position << nl
                << "which is in cell " << realCellI << nl
                << endl << abort(FatalError);
        }
        else
        {
            FatalErrorIn("gradInterpolationConstantTet<Type>::interpolate"
                        "("
                            "const vector&, "
                            "const label, "
                            "const label"
                        ") const")
                << "Failed to find tetrahedron containing point " << position << nl
                << "which is outside the domain.\n"
                << endl << abort(FatalError);
        }
    }
    return grad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
