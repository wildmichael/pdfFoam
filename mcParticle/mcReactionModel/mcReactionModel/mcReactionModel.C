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

#include "mcReactionModel.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcReactionModel, 0);
    defineRunTimeSelectionTable(mcReactionModel, mcReactionModel);

} // namespace Foamm

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcReactionModel::mcReactionModel
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcModel(cloud, db, subDictName)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mcReactionModel> Foam::mcReactionModel::New
(
    mcParticleCloud& cloud,
    const objectRegistry& db
)
{
    word reactionType(cloud.thermoDict().lookup("reactionModel"));
    word sd = reactionType + "ReactionModelCoeffs";

    mcReactionModelConstructorTable::iterator cstrIter =
        mcReactionModelConstructorTablePtr_->find(reactionType);

    if (cstrIter == mcReactionModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mcReactionModel::New(const objectRegistry&, const dictionary&)"
        )   << "Unknown mcReactionModel type " << reactionType << endl << endl
            << "Valid mcReactionModel types are :" << endl
            << mcReactionModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<mcReactionModel>
    (
        cstrIter()(cloud, db, sd)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::mcReactionModel::findIdx
(
    const word& nameKey,
    const word& defaultName
)
{
    word name = thermoDict().lookupOrDefault<word>(nameKey, defaultName);
    const wordList& sn = cloud().scalarNames();
    wordList::const_iterator f = sn.cbegin(), e = sn.cend();
    label idx = std::distance(f, std::find(f, e, name));
    if (idx == sn.size())
    {
        FatalErrorIn
        (
            "mcReactionModel::findIdx(const word&, const word&)"
        )
        << "Failed to find a scalar property names `" << name << "'\n"
        << exit(FatalError);
    }
    return idx;
}

// ************************************************************************* //
