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

#include "mcBurkeSchumannReactionModel.H"

#include "addToRunTimeSelectionTable.H"
#include "mcParticleCloud.H"

#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcBurkeSchumannReactionModel, 0);
    addNamedToRunTimeSelectionTable
    (
        mcReactionModel,
        mcBurkeSchumannReactionModel,
        mcReactionModel,
        BurkeSchumann
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcBurkeSchumannReactionModel::mcBurkeSchumannReactionModel
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& parentDict,
    const Foam::dictionary& mcBurkeSchumannReactionModelDict
)
:
    mcReactionModel(db, parentDict, mcBurkeSchumannReactionModelDict),
    Rfuel_(readScalar(lookup("Rfuel"))),
    Rox_  (readScalar(lookup("Rox"))),
    Rstoi_(readScalar(lookup("Rstoi"))),
    Tfuel_(readScalar(lookup("Tfuel"))),
    Tox_  (readScalar(lookup("Tox"))),
    Tstoi_(readScalar(lookup("Tstoi"))),
    zstoi_(readScalar(lookup("zstoi"))),
    zIdx_(-1),
    p_
    (
        db.lookupObject<volScalarField>
            (lookupOrDefault<word>("pName", "p", true))
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcBurkeSchumannReactionModel::correct(Foam::mcParticleCloud& cloud)
{
    forAllIter(mcParticleCloud, cloud, pIter)
    {
        correct(cloud, pIter());
    }
}


void Foam::mcBurkeSchumannReactionModel::correct
(
    Foam::mcParticleCloud& cloud,
    Foam::mcParticle& p
)
{
    setTIdx(cloud);
    if (zIdx_ < 0)
    {
        word zName = lookupOrDefault<word>("zName", "z", true);
        const wordList& sn = cloud.scalarNames();
        wordList::const_iterator f = sn.cbegin(), e = sn.cend();
        zIdx_ = std::distance(f, std::find(f, e, zName));
        if (zIdx_ == sn.size())
        {
            FatalErrorIn
            (
                "mcBurkeSchumannReactionModel::correct("
                "mcParticleCloud&, mcParticle&)"
            )
            << "Failed to find a scalar property names `" << zName << "'\n"
            << exit(FatalError);
        }
    }
    const scalar& z = p.Phi()[zIdx_];
    static const scalar RTox = Rox_*Tox_;
    static const scalar RTfuel = Rfuel_*Tfuel_;
    static const scalar RTstoi = Rstoi_*Tstoi_;
    scalar RT = z<zstoi_
        ? (RTstoi - RTox)/zstoi_*z + RTox
        : (RTfuel - RTstoi)/(1 - zstoi_)*(z - zstoi_) + RTstoi;
    scalar T = z<zstoi_
        ? (Tstoi_ - Tox_)/zstoi_*z + Tox_
        : (Tfuel_ - Tstoi_)/(1 - zstoi_)*(z - zstoi_) + Tstoi_;
    p.rho() = p_[p.cell()]/(RT);
    p.Phi()[TIdx_] = T;
}

// ************************************************************************* //
