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

#include "tetFacePointCellDecomposition.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Tetrahedron>
Foam::tetFacePointCellDecomposition<Tetrahedron>
::tetFacePointCellDecomposition(const polyMesh& pMesh)
:
    pMesh_(pMesh),
    tets_(decompose(pMesh_)().xfer()),
    cellTets_(pMesh.cells().size()),
    tetFace_(tets_.size()),
    tetPoints_(tets_.size())
{
    const cellList& cells = this->pMesh_.cells();
    const faceList& faces = this->pMesh_.faces();
    label tetI = 0;
    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];
        DynamicList<label> ctets;
        forAll(c, faceI)
        {
            label gFaceI = c[faceI];
            const face& f = faces[gFaceI];
            forAll(f, i)
            {
                ctets.append(tetI);
                tetFace_[tetI] = gFaceI;
                tetPoints_[tetI].first() = i;
                tetPoints_[tetI].second() = f.rcIndex(i);
                if (cellI != this->pMesh_.faceOwner()[gFaceI])
                {
                    tetPoints_[tetI] = tetPoints_[tetI].reversePair();
                }
                ++tetI;
            }
        }
        cellTets_[cellI].transfer(ctets);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Tetrahedron>
Foam::autoPtr<Foam::List<Tetrahedron> >
Foam::tetFacePointCellDecomposition<Tetrahedron>::decompose
(
    const polyMesh& pMesh
)
{
    typedef List<Tetrahedron> tetList;
    const cellList& cells = pMesh.cells();
    const faceList& pMeshFaces = pMesh.faces();
    label nTets = 0;
    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];
        forAll(c, cellFaceI)
        {
            nTets += pMeshFaces[c[cellFaceI]].size();
        }
    }
    autoPtr<List<Tetrahedron> > tetsPtr(new List<Tetrahedron>(nTets));
    List<Tetrahedron>& tets = tetsPtr();
    label tetI = 0;
    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];
        forAll(c, cellFaceI)
        {
            const face& f = pMeshFaces[c[cellFaceI]];
            forAll(f, facePointI)
            {
                if (tetI==nTets)
                {
                    FatalErrorIn("tetFacePointCellDecomposition::decompose()")
                        << "Somehow miscounted tetrahedra in mesh.\n"
                        << endl << abort(FatalError);
                }
                tets[tetI++] =
                    tetFromFacePointsCell
                    (
                        pMesh,
                        cellI,
                        cellFaceI,
                        facePointI
                    )();
            }
        }
    }
    return tetsPtr;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Tetrahedron>
Foam::label Foam::tetFacePointCellDecomposition<Tetrahedron>::find
(
    const point& pt,
    label cellHint
) const
{
    if (cellHint > -1)
    {
        const labelList& ctets = cellTets_[cellHint];
        forAll(ctets, cteti)
        {
            label teti = ctets[cteti];
            if (tets_[teti].inside(pt))
            {
                return teti;
            }
        }
    }
    forAll(tets_, teti)
    {
        if (tets_[teti].inside(pt))
        {
            return teti;
        }
    }
    return -1;
}


// ************************************************************************* //
