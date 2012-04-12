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

#include "mcNaiveFlameletModel.H"

#include "addToRunTimeSelectionTable.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcNaiveFlameletModel, 0);
    addNamedToRunTimeSelectionTable
    (
        mcReactionModel,
        mcNaiveFlameletModel,
        mcReactionModel,
        naiveFlamelet
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcNaiveFlameletModel::mcNaiveFlameletModel
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& parentDict,
    const Foam::dictionary& mcNaiveFlameletModelDict
)
:
    mcReactionModel(db, parentDict, mcNaiveFlameletModelDict),
    zName_(lookupOrDefault<word>("zName", "z", true)),
    zIdx_(-1)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcNaiveFlameletModel::correct(Foam::mcParticleCloud& cloud)
{
    forAllIter(mcParticleCloud, cloud, pIter)
    {
        correct(cloud, pIter());
    }
}


void Foam::mcNaiveFlameletModel::correct
(
    Foam::mcParticleCloud& cloud,
    Foam::mcParticle& p
)
{
    if (zIdx_ < 0)
    {
        const wordList& names = cloud.scalarNames();
        wordList::const_iterator zIter = find(names.begin(), names.end(),
                                              zName_);
        if (zIter == names.end())
        {
            Foam::FatalErrorIn
            (
                "mcNaiveFlameletModel::correct(mcParticleCloud&, mcParticle&)"
            )
                << "No such particle property: " << zName_ << endl
                << exit(FatalError);
        }
        zIdx_ = distance(names.begin(), zIter);
    }
    const scalar& z = p.Phi()[zIdx_];
    p.rho() = (1. - 3.2*z*(1.-z));
}

// ************************************************************************* //
