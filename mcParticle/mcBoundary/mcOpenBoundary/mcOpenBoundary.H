/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::mcOpenBoundary

Description
    Base-class for open Monte Carlo particle boundaries.

    If the @c reflecting switch is enabled (the default), consistent reflecting
    is enabled (see D.W. Meyer and P. Jenny, JCP 226 1859-1873,
    http://dx.doi.org/10.1016/j.jcp.2007.06.014).

    The reflection at the inflow boundary is insofar inaccurate as the
    wall-parallel motion outside the domain is neglected.

SourceFiles
    mcOpenBoundary.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcOpenBoundary_H
#define mcOpenBoundary_H

#include "mcBoundary.H"

#include "vectorList.H"
#include "PtrMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class mcOpenBoundary Declaration
\*---------------------------------------------------------------------------*/

class mcOpenBoundary :
    public mcBoundary
{
    // Private Types

        typedef PtrMap<DynamicList<label> > dynamicLabelListPtrMap;

    // Private data

        //- Is this boundary reflecting?
        Switch reflecting_;

        //- Mapping of cells adjacent to this patch to a list of face indices
        dynamicLabelListPtrMap cellFaces_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcOpenBoundary(const mcOpenBoundary&);
        void operator=(const mcOpenBoundary&);

protected:

    // Protected data

        //- Patch normals
        vectorField n_;

        //- Patch-normal velocity field
        scalarField Un_;

public:

    //- Runtime type information
    TypeName("mcOpenBoundary");


    // Constructors

        //- Construct from components
        mcOpenBoundary
        (
            mcParticleCloud& cloud,
            label patchID,
            const dictionary& dict
        );

    // Member Functions

        //- Whether boundary is reflecting
        const Switch& reflecting() const { return reflecting_; }

        //- Whether boundary is reflecting
        Switch& reflecting() { return reflecting_; }

        //- Delete particles at the outlet boundary
        virtual void correct(bool afterMove);

        //- Deal with a particle hitting the patch
        virtual void hitPatch
        (
            mcParticle& p,
            mcParticle::trackData& td
        );

        //- Deal with a particle hitting the patch without trackData
        virtual void hitPatch
        (
            mcParticle& p,
            int&
        );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
