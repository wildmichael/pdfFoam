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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::mcParticleCloud& Foam::mcBoundary::cloud() const
{
    return cloud_;
}


Foam::mcParticleCloud& Foam::mcBoundary::cloud()
{
    return cloud_;
}


Foam::label Foam::mcBoundary::patchID() const
{
    return patchID_;
}


const Foam::polyPatch& Foam::mcBoundary::patch() const
{
    return patch_;
}


const Foam::scalarField& Foam::mcBoundary::scalarInFlux() const
{
    return scalarInFlux_;
}


Foam::scalarField& Foam::mcBoundary::scalarInFlux()
{
    return scalarInFlux_;
}


const Foam::scalarField& Foam::mcBoundary::scalarOutFlux() const
{
    return scalarOutFlux_;
}


Foam::scalarField& Foam::mcBoundary::scalarOutFlux()
{
    return scalarOutFlux_;
}

// ************************************************************************* //
