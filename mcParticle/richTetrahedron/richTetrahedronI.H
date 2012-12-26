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

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Foam::richTetrahedron<Point, PointRef>::richTetrahedron()
:
    tet_
    (
        new tetrahedron<Point, PointRef>
        (
            point::zero,
            point::zero,
            point::zero,
            point::zero
        )
    ),
    gradNi_(tetrahedron<Point, PointRef>::nVertices)
{
    update();
}

template<class Point, class PointRef>
inline Foam::richTetrahedron<Point, PointRef>::richTetrahedron
(
    const Point& a,
    const Point& b,
    const Point& c,
    const Point& d
)
:
    tet_(new tetrahedron<Point, PointRef>(a, b, c, d)),
    gradNi_(tetrahedron<Point, PointRef>::nVertices)
{
    update();
}


template<class Point, class PointRef>
inline Foam::richTetrahedron<Point, PointRef>::richTetrahedron(Istream& is)
    : tet_(new tetrahedron<Point, PointRef>(is)),
    gradNi_(tetrahedron<Point, PointRef>::nVertices)
{
    update();
}


template<class Point, class PointRef>
inline Foam::richTetrahedron<Point, PointRef>::richTetrahedron
(
    const richTetrahedron<Point, PointRef>& other
) :
    tet_(new tetrahedron<Point, PointRef>(other.tet_())),
    Sa_(other.Sa_),
    Sb_(other.Sb_),
    Sc_(other.Sc_),
    Sd_(other.Sd_),
    mag_(other.mag_),
    gradNi_(other.gradNi_)
{}


template<class Point, class PointRef>
inline void Foam::richTetrahedron<Point, PointRef>::operator=
(
    const richTetrahedron<Point, PointRef>& other
)
{
    tet_.reset(new tetrahedron<Point, PointRef>(other.tet_()));
    Sa_ = other.Sa_;
    Sb_ = other.Sb_;
    Sc_ = other.Sc_;
    Sd_ = other.Sd_;
    mag_ = other.mag_;
    gradNi_ = other.gradNi_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline void Foam::richTetrahedron<Point, PointRef>::update()
{
    Sa_  = tet_().Sa();
    Sb_  = tet_().Sb();
    Sc_  = tet_().Sc();
    Sd_  = tet_().Sd();
    mag_ = tet_().mag();
    if (mag_ > 0)
    {
        gradNi_[0] = Sa_;
        gradNi_[1] = Sb_;
        gradNi_[2] = Sc_;
        gradNi_[3] = Sd_;
        gradNi_ /= -3.*mag_;
    }
    else
    {
        gradNi_ = vector::zero;
    }
#ifdef FULLDEBUG
    if (mag_ < 0)
    {
        WarningIn("richTetrahedron::update()")
            << "Tetrahedron has negative volume: " << mag_ << nl
            << endl;
    }
#endif
}


template<class Point, class PointRef>
inline Foam::vector Foam::richTetrahedron<Point, PointRef>::Sa() const
{
    return Sa_;
}


template<class Point, class PointRef>
inline Foam::vector Foam::richTetrahedron<Point, PointRef>::Sb() const
{
    return Sb_;
}


template<class Point, class PointRef>
inline Foam::vector Foam::richTetrahedron<Point, PointRef>::Sc() const
{
    return Sc_;
}


template<class Point, class PointRef>
inline Foam::vector Foam::richTetrahedron<Point, PointRef>::Sd() const
{
    return Sd_;
}


template<class Point, class PointRef>
inline Foam::scalar Foam::richTetrahedron<Point, PointRef>::mag() const
{
    return mag_;
}


template<class Point, class PointRef>
inline const Foam::vectorField&
Foam::richTetrahedron<Point, PointRef>::gradNi() const
{
    return gradNi_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Foam::Istream&
Foam::operator>>(Istream& is, richTetrahedron<Point, PointRef>& t)
{
    is >> t.tet_();
    t.update();
    return is;
}


template<class Point, class PointRef>
inline Foam::Ostream&
Foam::operator<<(Ostream& os, const richTetrahedron<Point, PointRef>& t)
{
    return (os << t.tet_());
}


template<class Point, class PointRef>
inline bool
Foam::operator!=
(
    const richTetrahedron<Point, PointRef>& t1,
    const richTetrahedron<Point, PointRef>& t2
)
{
    return t1.a() != t2.a()
        || t1.b() != t2.b()
        || t1.c() != t2.c()
        || t1.d() != t2.d();
}


template<class Point, class PointRef>
inline bool
Foam::operator==
(
    const richTetrahedron<Point, PointRef>& t1,
    const richTetrahedron<Point, PointRef>& t2
)
{
    return !operator!=(t1, t2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
