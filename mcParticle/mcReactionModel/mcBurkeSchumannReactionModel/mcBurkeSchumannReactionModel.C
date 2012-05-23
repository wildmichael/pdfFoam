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
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcReactionModel(cloud, db, subDictName),
    Rfuel_(readScalar(thermoDict().lookup("Rfuel"))),
    Rox_  (readScalar(thermoDict().lookup("Rox"))),
    Rstoi_(readScalar(thermoDict().lookup("Rstoi"))),
    Tfuel_(readScalar(thermoDict().lookup("Tfuel"))),
    Tox_  (readScalar(thermoDict().lookup("Tox"))),
    Tstoi_(readScalar(thermoDict().lookup("Tstoi"))),
    zstoi_(readScalar(thermoDict().lookup("zstoi"))),
    zIdx_(findIdx("zName", "z")),
    TIdx_(findIdx("TName", "T")),
    p_
    (
        db.lookupObject<volScalarField>
            (thermoDict().lookupOrDefault<word>("pName", "p"))
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcBurkeSchumannReactionModel::correct(mcParticle& p)
{
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
