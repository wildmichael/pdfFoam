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
    Foam::gradInterpolationConstantTet

Description
    Decomposes the FV cell into tetrahedrons consisting of face centres, cell
    centres and points which build linear finite elements. The interpolated
    gradient is then piecewise constant within the tetrahedrons.

SourceFiles
    gradInterpolationConstantTet.C
    makeGradInterpolationConstantTet.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef gradInterpolationConstantTet_H
#define gradInterpolationConstantTet_H

#include "gradInterpolation.H"
#include "richTetPointRef.H"
#include "tetFacePointCellDecomposition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class gradInterpolationConstantTet Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gradInterpolationConstantTet
:
    public gradInterpolation<Type>
{
    // Private data

        //- Interpolated volfield
        const GeometricField<Type, pointPatchField, pointMesh> psip_;

        //- Linearly interpolated volfield
        const GeometricField<Type, fvsPatchField, surfaceMesh> psis_;

        //- The mesh decomposition
        tmp<tetFacePointCellDecomposition<richTetPointRef> > tetDecomp_;


public:

    //- Runtime type information
    TypeName("constantTet");


    // Constructors

        //- Construct from components
        gradInterpolationConstantTet
        (
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        //- Construct given a tet decomposition
        gradInterpolationConstantTet
        (
            const tetFacePointCellDecomposition<richTetPointRef>& tetDecomp,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );


    // Member Functions

        //- Interpolate field to the given point in the given cell
        // N.B: outerProduct<vector, Type>::type computes the correct return
        // type at compile time (gradient of a scalar is a vector, gradient of
        // a vector is a tensor, etc.).
        typename outerProduct<vector, Type>::type interpolate
        (
            const vector& position,
            const label cellI,
            const label faceI = -1
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gradInterpolationConstantTet.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
