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
#include "tetrahedronTools.H"

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
    tets_
    (
        tetrahedronTools::decomposePointFaceCell<richTetPointRef>
        (
            psi.mesh()
        ).xfer()
    ),
    cellTets_(this->pMesh_.cells().size())
{
    const cellList& cells = this->pMesh_.cells();
    label tetI = 0;
    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];
        label nTets = 0;
        forAll(c, faceI)
        {
            nTets += this->pMeshFaces_[c[faceI]].size();
        }
        cellTets_[cellI].reset
        (
            new SubList<richTetPointRef>
            (
                tets_,
                nTets,
                tetI
            )
        );
        tetI += nTets;
    }
}


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
    const cell& c = this->pMesh_.cells()[cellI];
    typename Foam::outerProduct<Foam::vector, Type>::type grad;
    label tetI =
        tetrahedronTools::findTetrahedron(cellTets_[cellI](), position);
    if (tetI > -1)
    {
        // find the face this tet "stands" on
        label n = 0;
        label cellFaceI = 0;
        label facePointI = -1;
        for (; cellFaceI != c.size(); ++cellFaceI)
        {
            const face& f = this->pMeshFaces_[c[cellFaceI]];
            label s = f.size();
            if (tetI < n+s)
            {
                facePointI = tetI - n;
                break;
            }
            n += s;
        }
        if (facePointI < 0)
        {
            FatalErrorIn("gradInterpolationConstantTet<Type>::interpolate"
                "(const vector&, const label, const label)")
                << "Failed to find face and point belonging to tetrahedron\n"
                << endl << abort(FatalError);
        }
        label gFaceI = c[cellFaceI];
        const face& f = this->pMeshFaces_[gFaceI];
#ifdef FULLDEBUG
        if (facePointI >= f.size())
        {
            FatalErrorIn("gradInterpolationConstantTet<Type>::interpolate"
                "(const vector&, const label, const label)")
                << "Computed facePointI is larger than face size\n"
                << endl << abort(FatalError);
        }
#endif
        label ptBI = f[facePointI];
        label ptCI = f[f.rcIndex(facePointI)];
        if (cellI != pMesh.faceOwner()[faceI])
        {
            Swap(ptBI, ptCI);
        }
        const vectorField& gradNi =
            cellTets_[cellI]()[tetI].gradNi();
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
        FatalErrorIn("gradInterpolationConstantTet<Type>::interpolate"
                    "("
                        "const vector&, "
                        "const label, "
                        "const label"
                    ") const")
            << "Failed to find tetrahedron containing point " << position << nl
            << "The cell " << cellI << ": "
            << c.points(this->pMeshFaces_, this->pMeshPoints_) << nl
            << "The tetrahedra are: " << cellTets_[cellI]() << nl
            << endl << abort(FatalError);
    }
    return grad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
