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

#include "mcLocalTimeStepping.H"

#include "mcParticleCloud.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcLocalTimeStepping, 0);
    defineRunTimeSelectionTable(mcLocalTimeStepping, mcLocalTimeStepping);

} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcLocalTimeStepping::mcLocalTimeStepping
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& parentDict,
    const Foam::dictionary& mcLocalTimeSteppingDict
)
:
    mcModel(db, parentDict, mcLocalTimeSteppingDict),
    active_(true)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mcLocalTimeStepping> Foam::mcLocalTimeStepping::New
(
    const Foam::objectRegistry& db,
    const Foam::dictionary& dict
)
{
    word name(dict.lookup("localTimeStepping"));

    // If set to "false", "no", "n", "off" or "none", disable
    Switch::switchType enable = Switch::asEnum(name, true);
    if (enable != Switch::INVALID && !Switch::asBool(enable))
    {
        autoPtr<mcLocalTimeStepping> result
        (
            new mcLocalTimeStepping(db, dict, dictionary::null)
        );
        result().active_ = false;
        return result;
    }

    mcLocalTimeSteppingConstructorTable::iterator cstrIter =
        mcLocalTimeSteppingConstructorTablePtr_->find(name);

    if (cstrIter == mcLocalTimeSteppingConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mcLocalTimeStepping::New(const fvMesh&, const dictionary&)"
        )   << "Unknown mcLocalTimeStepping type " << name << nl << nl
            << "Valid mcLocalTimeStepping types are :" << nl
            << mcLocalTimeSteppingConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<mcLocalTimeStepping>
    (
        cstrIter()
        (
            db,
            dict,
            dict.subOrEmptyDict(name+"LocalTimeSteppingCoeffs")
        )
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcLocalTimeStepping::setupInternals(const mcParticleCloud& cloud)
{}


void Foam::mcLocalTimeStepping::correct(mcParticleCloud& cloud)
{
    setupInternals(cloud);
    forAllIter(mcParticleCloud, cloud, pIter)
    {
        this->correct(cloud, pIter(), false);
    }
}


void Foam::mcLocalTimeStepping::correct
(
    mcParticleCloud&,
    mcParticle&,
    bool prepare
)
{}

// ************************************************************************* //
