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
    Foam::gradInterpolation

Description
    Abstract base class for gradInterpolation. These interpolation classes are
    intended to take values at the cell centres and interpolate the gradient of
    these values to the given location.

SourceFiles
    gradInterpolation.C
    gradInterpolationNew.C
    gradInterpolations.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef gradInterpolation_H
#define gradInterpolation_H

#include "faceList.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class gradInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gradInterpolation
{

protected:

    // Protected data

        const GeometricField<Type, fvPatchField, volMesh>& psi_;

        const polyMesh& pMesh_;
        const vectorField& pMeshPoints_;
        const faceList& pMeshFaces_;
        const vectorField& pMeshFaceCentres_;
        const vectorField& pMeshFaceAreas_;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            gradInterpolation,
            dictionary,
            (
                const GeometricField<Type, fvPatchField, volMesh>& psi
            ),
            (psi)
        );


    // Selectors

        //- Return a reference to the specified gradInterpolation scheme
        static autoPtr<gradInterpolation<Type> > New
        (
            const word& gradInterpolationType,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        //- Return a reference to the selected gradInterpolation scheme
        static autoPtr<gradInterpolation<Type> > New
        (
            const dictionary& interpolationSchemes,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );


    // Constructors

        //- Construct from components
        gradInterpolation
        (
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );


    //- Destructor
    virtual ~gradInterpolation()
    {}


    // Member Functions

        //- Return the field to be interpolated
        const GeometricField<Type, fvPatchField, volMesh>& psi() const
        {
            return psi_;
        }

        //- Interpolate gradient of field to the given point in the given cell
        virtual typename outerProduct<vector, Type>::type interpolate
        (
            const vector& position,
            const label cellI,
            const label faceI = -1
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeGradInterpolationType(SS, Type)                                    \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
gradInterpolation<Type>::adddictionaryConstructorToTable<SS<Type> >            \
    add##SS##Type##ConstructorToTable_;


#define makeGradInterpolation(SS)                                              \
                                                                               \
makeGradInterpolationType(SS, scalar)                                          \
makeGradInterpolationType(SS, vector)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gradInterpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
