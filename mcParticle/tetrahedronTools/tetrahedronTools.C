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

#include "tetrahedronTools.H"

// * * * * * * * * * * * * * Local Helper Functions  * * * * * * * * * * * * //

namespace // anonymous
{


} // anonymous namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Tetrahedron>
inline Foam::autoPtr<Tetrahedron> Foam::tetrahedronTools::tetFromPointsFaceCell
(
    const polyMesh& pMesh,
    label cellI,
    label cellFaceI,
    label facePointI
)
{
    const pointField& pMeshPoints = pMesh.points();
    const cell& c = pMesh.cells()[cellI];
    label faceI = c[cellFaceI];
    const face& f = pMesh.faces()[faceI];
    bool own = cellI == pMesh.faceOwner()[faceI];

    const point& ptA = pMesh.faceCentres()[faceI];
    const point& ptB = pMeshPoints[f[facePointI]];
    const point& ptC = pMeshPoints[f[f.rcIndex(facePointI)]];
    const point& ptD = pMesh.cellCentres()[cellI];

    autoPtr<Tetrahedron> tetPtr;
    if (own)
    {
        tetPtr.reset(new Tetrahedron(ptA, ptB, ptC, ptD));
    }
    else
    {
        tetPtr.reset(new Tetrahedron(ptA, ptC, ptB, ptD));
    }

    return tetPtr;
}


template<class Tetrahedron>
inline Foam::List<Tetrahedron>
Foam::tetrahedronTools::decomposePointFaceCell
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
    List<Tetrahedron> tets(nTets);
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
                    FatalErrorIn("tetrahedronTools::decomposePointFaceCell"
                        "(const polyMesh&)")
                        << "Somehow miscounted tetrahedra in mesh.\n"
                        << endl << abort(FatalError);
                }
                tets[tetI++] =
                    tetFromPointsFaceCell<Tetrahedron>
                    (
                        pMesh,
                        cellI,
                        cellFaceI,
                        facePointI
                    )();
            }
        }
    }
    return tets;
}


template<class Tetrahedron>
inline Foam::label Foam::tetrahedronTools::findTetrahedron
(
    const UList<Tetrahedron>& tets,
    const point& pt
)
{
    forAll(tets, tetI)
    {
        if (tets[tetI].inside(pt))
        {
            return tetI;
        }
    }
    return -1;
}

// ************************************************************************* //
