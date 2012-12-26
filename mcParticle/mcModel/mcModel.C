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

#include "mcModel.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcModel, 0);

} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcModel::mcModel
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    cloud_(cloud),
    db_(db),
    subDictName_(subDictName),
    thermoDict_(cloud.thermoDict().subOrEmptyDict(subDictName_))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary Foam::mcModel::solutionDict() const
{
    const mcSolution& s = cloud().solutionDict();
    if (s.found(subDictName_) && !s.solutionDict().found(subDictName_))
    {
        // If the user used "select", but did not provide a model-specific
        // sub-dictionary, revert to the top-level one.
        return s.subOrEmptyDict(subDictName_);
    }
    return s.solutionDict()
        .subOrEmptyDict(subDictName_);
}


void Foam::mcModel::updateInternals()
{}


void Foam::mcModel::correct()
{
    updateInternals();
    forAllIter(mcParticleCloud, cloud_, pIter)
    {
        correct(pIter());
    }
}


void Foam::mcModel::Co(mcParticle&) const
{}

// ************************************************************************* //
