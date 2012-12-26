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

#include "mcEmptyBoundary.H"

#include "addToRunTimeSelectionTable.H"
#include "mcParticle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcEmptyBoundary, 0);
    addNamedToRunTimeSelectionTable
    (
        mcBoundary,
        mcEmptyBoundary,
        mcBoundary,
        empty
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcEmptyBoundary::mcEmptyBoundary
(
    mcParticleCloud& cloud,
    label patchID,
    const dictionary& dict
)
:
    mcBoundary(cloud, patchID, dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#if FOAM_HEX_VERSION < 0x200

void Foam::mcEmptyBoundary::hitPatch
(
    mcParticle& p,
    mcParticle::trackData& td
)
{
    FatalErrorIn("mcEmptyBoundary::hitPatch"
        "(mcParticle&, mcParticle::trackData&)")
        << "A mcParticle must never hit an empty boundary!" << endl
        << Foam::exit(FatalError);
}


void Foam::mcEmptyBoundary::hitPatch
(
    mcParticle& p,
    int&
)
{
    FatalErrorIn("mcEmptyBoundary::hitPatch"
        "(mcParticle&, int&)")
        << "A mcParticle must never hit an empty boundary!" << endl
        << Foam::exit(FatalError);
}

#else

void Foam::mcEmptyBoundary::hitPatch
(
    mcParticle& p,
    mcParticle::trackData& td,
    const label patchI,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    FatalErrorIn("mcEmptyBoundary::hitPatch"
        "(mcParticle&, mcParticle::trackData&, const label, const scalar, "
        "const tetIndices&)")
        << "A mcParticle must never hit an empty boundary!" << endl
        << Foam::exit(FatalError);
}

#endif

// ************************************************************************* //
