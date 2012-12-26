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

Class
    Foam::richTetrahedron

Description
    A tetrahedron built from two adjacent points, the face centre and the cell
    centre. In contrast to Foam::tetrahedron, this class stores the face
    normals and adds additional functions.

SourceFiles
    richTetrahedronI.H
    richTetrahedron.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef richTetrahedron_H
#define richTetrahedron_H

#include "tetrahedron.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators

template<class Point, class PointRef> class richTetrahedron;

template<class Point, class PointRef>
inline Istream& operator>>
(
    Istream&,
    richTetrahedron<Point, PointRef>&
);

template<class Point, class PointRef>
inline Ostream& operator<<
(
    Ostream&,
    const richTetrahedron<Point, PointRef>&
);

template<class Point, class PointRef>
inline bool operator!=
(
    const richTetrahedron<Point, PointRef>&,
    const richTetrahedron<Point, PointRef>&
);

template<class Point, class PointRef>
inline bool operator==
(
    const richTetrahedron<Point, PointRef>&,
    const richTetrahedron<Point, PointRef>&
);


/*---------------------------------------------------------------------------*\
                  Class richTetrahedron Declaration
\*---------------------------------------------------------------------------*/

template<class Point, class PointRef>
class richTetrahedron
{
    // Private data

        autoPtr<tetrahedron<Point, PointRef> > tet_;

        //- Surface vector of face opposite of point a_
        vector Sa_;

        //- Surface vector of face opposite of point b_
        vector Sb_;

        //- Surface vector of face opposite of point c_
        vector Sc_;

        //- Surface vector of face opposite of point d_
        vector Sd_;

        //- Volume
        scalar mag_;

        //- Gradients of shape functions
        vectorField gradNi_;

    // Private Member Functions

        //- Update surface normals, volume and shape function gradients
        inline void update();

public:


    // Constructors

        //- Construct null
        inline richTetrahedron();

        //- Construct from points
        inline richTetrahedron
        (
            const Point& a,
            const Point& b,
            const Point& c,
            const Point& d
        );

        //- Construct from Istream
        explicit richTetrahedron(Istream&);

        //- Construct as a copy
        explicit richTetrahedron(const richTetrahedron& other);

        //- Construct and return a clone
        autoPtr<richTetrahedron<Point, PointRef> > clone() const
        {
            return autoPtr<richTetrahedron<Point, PointRef> >
            (
                new richTetrahedron<Point, PointRef>(*this)
            );
        }

        //- Assignment operator
        void operator=(const richTetrahedron& other);


    // Member Functions

        void swap(richTetrahedron& other);

        // Access

            const Point& a() const {return tet_().a();}

            const Point& b() const {return tet_().b();}

            const Point& c() const {return tet_().c();}

            const Point& d() const {return tet_().d();}


        // Properties

            //- Return face normal
            inline vector Sa() const;

            inline vector Sb() const;

            inline vector Sc() const;

            inline vector Sd() const;

            Point centre() const {return tet_().centre();}

            //- Return volume
            inline scalar mag() const;

            Point circumCentre() const {return tet_().circumCentre();}

            scalar circumRadius() const {return tet_().circumRadius();}

            pointHit containmentSphere(const scalar tol) const
            {return tet_().containmentSphere(tol);}

            //- Gradients of shape functions
            const vectorField& gradNi() const;

            //- Return true if point is inside tetrahedron
            bool inside(const point& pt) const;

    // Friend operators

        friend Istream& operator>> <Point, PointRef>
        (
            Istream&,
            richTetrahedron&
        );

        friend Ostream& operator<< <Point, PointRef>
        (
            Ostream&,
            const richTetrahedron&
        );

        friend bool operator!= <Point, PointRef>
        (
            const richTetrahedron& t1,
            const richTetrahedron& t2
        );

        friend bool operator== <Point, PointRef>
        (
            const richTetrahedron& t1,
            const richTetrahedron& t2
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "richTetrahedronI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "richTetrahedron.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
