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

#include "mcVelocityModel.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcVelocityModel, 0);
    defineRunTimeSelectionTable(mcVelocityModel, mcVelocityModel);

} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcVelocityModel::mcVelocityModel
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& parentDict,
    const Foam::dictionary& dict,
    Foam::mcParticleCloud& cloud
)
:
    mcModel(db, parentDict, dict),
    cloud_(cloud)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mcVelocityModel> Foam::mcVelocityModel::New
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& dict,
    Foam::mcParticleCloud& cloud
)
{
    word velocityModelType(dict.lookup("velocityModel"));

    mcVelocityModelConstructorTable::iterator cstrIter =
        mcVelocityModelConstructorTablePtr_->find(velocityModelType);

    if (cstrIter == mcVelocityModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mcVelocityModel::New(const objectRegistry&, const dictionary&)"
        )   << "Unknown mcVelocityModel type " << velocityModelType
            << endl << endl
            << "Valid mcVelocityModel types are :" << endl
            << mcVelocityModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<mcVelocityModel>
    (
        cstrIter()
        (
            db,
            dict,
            dict.subOrEmptyDict(velocityModelType+"VelocityModelCoeffs"),
            cloud
        )
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcVelocityModel::setupInternals()
{}


void Foam::mcVelocityModel::correct()
{
    setupInternals();
    forAllIter(mcParticleCloud, cloud_, pIter)
    {
        correct(pIter(), false);
    }
}

// ************************************************************************* //
