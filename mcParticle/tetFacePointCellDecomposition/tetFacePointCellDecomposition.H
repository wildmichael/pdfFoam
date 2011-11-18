/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::tetFacePointCellDecomposition

Description
    Decomposes a polyMesh into tetrahedra consisting of the face centre two
    points and the cell centre.

SourceFiles
    tetFacePointCellDecompositionI.H
    tetFacePointCellDecomposition.C
    tetFacePointCellDecompositionIO.C

\*---------------------------------------------------------------------------*/

#ifndef tetFacePointCellDecomposition_H
#define tetFacePointCellDecomposition_H

#include "autoPtr.H"
#include "labelList.H"
#include "labelPair.H"
#include "refCount.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;


/*---------------------------------------------------------------------------*\
                Class tetFacePointCellDecomposition Declaration
\*---------------------------------------------------------------------------*/

template<class Tetrahedron>
class tetFacePointCellDecomposition
:
    public refCount
{
    // Private data

        //- The mesh
        const polyMesh& pMesh_;

        //- The tetrahedra
        List<Tetrahedron> tets_;

        //- The tetrahedra by cell
        labelListList cellTets_;

        //- The face a tetrahedron belongs to
        labelList tetFace_;

        //- The two face points in a tetrahedron
        List<labelPair> tetPoints_;


    // Private Static Member Functions

        //- Return new tetrahedron formed by two points, face and cell centres
        static inline autoPtr<Tetrahedron> tetFromFacePointsCell
        (
            const polyMesh& pMesh, //!< The mesh
            label cellI,     //!< Global cell index containing the tetrahedron
            label cellFaceI, //!< Index of the cell-face containing the base
            label facePointI //!< Index of the first point in the face
        );

        //- Decompose mesh into tetrahedra of two points, face and cell centres
        // Returns a list containing the tetrahedra. The ordering is first by
        // cell, then by faces and finally by points.
        static autoPtr<List<Tetrahedron> > decompose(const polyMesh& pMesh);


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetFacePointCellDecomposition(const tetFacePointCellDecomposition<Tetrahedron>&);

        //- Disallow default bitwise assignment
        void operator=(const tetFacePointCellDecomposition<Tetrahedron>&);


public:


    // Constructors

        //- Construct from components
        tetFacePointCellDecomposition(const polyMesh& pMesh);


    // Member Functions

        //- The list of tetrahedra
        inline const UList<Tetrahedron>& tetrahedra() const;

        //- The list of tetrahedra in every cell
        inline const labelListList& cellTetrahedra() const;

        //- The face a tetrahedron stands on
        inline const labelList& tetrahedronFace() const;

        //- The two face points in a tetrahedron
        inline const List<labelPair>& tetrahedronPoints() const;

        //- Find the tetrahedron containing the point
        // \returns The index of the tetrahedron, -1 if not found
        label find
        (
            const point& pt,     //!< Point to find
            label cellHint = -1  //!< Hint in which cell to look first
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetFacePointCellDecompositionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tetFacePointCellDecomposition.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
