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

Class

Description
    Utility functions for tetrahedrons

SourceFiles
    tetrahedronTools.C

\*---------------------------------------------------------------------------*/

#ifndef tetrahedronTools_H
#define tetrahedronTools_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class tetrahedronTools
{
public:

    //- Return new tetrahedron formed by two points, face and cell centres
    template<class Tetrahedron>
    static inline autoPtr<Tetrahedron> tetFromPointsFaceCell
    (
        const polyMesh& pMesh, //!< The mesh
        label cellI, //!< Global cell index containing the tetrahedron
        label cellFaceI, //!< Index of the cell-face containing the base
        label facePointI //!< Index of the first point in the face
    );


    //- Decompose mesh into tetrahedra of two points, face and cell centres
    // Returns a list containing the tetrahedra. The ordering is first by
    // cell, then by faces and finally by points.
    template<class Tetrahedron>
    static inline List<Tetrahedron> decomposePointFaceCell
    (
        const polyMesh& pMesh
    );


    //- Find the tetrahedron containing the point
    // \returns The index of the tetrahedron, -1 if not found
    template<class Tetrahedron>
    static inline label findTetrahedron
    (
        const UList<Tetrahedron>& tets,
        const point& pt
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tetrahedronTools.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
