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

#include <labelList.H>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- @todo This is a hack to fix the implementation in OpenFOAM < 2.0.x
// @see http://openfoam.org/mantisbt/view.php?id=330
template<class T>
void uniqueOrder_FIX
(
    const Foam::UList<T>& lst,
    Foam::labelList& order
)
{
    using namespace Foam;
    sortedOrder(lst, order);

    if (order.size() > 1)
    {
        label n = 0;
        for (label i = 0; i < order.size() - 1; ++i)
        {
            if (lst[order[i]] != lst[order[i+1]])
            {
                order[n++] = order[i];
            }
        }
        // DONT FORGET THE LAST ELEMENT!
        order[n++] = order[order.size()-1];
        order.setSize(n);
    }
}

// ************************************************************************* //
