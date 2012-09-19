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

#include "mcSolution.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcSolution::mcSolution(const objectRegistry& obr)
:
    IOdictionary
    (
        IOobject
        (
            "mcSolution",
            obr.time().system(),
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    averagingTime_("averagingTime", 0.01*obr.time().endTime()),
    defaultRelaxationTime_("relaxationTime", dimTime, GREAT),
    relaxationTimes_(ITstream("relaxationTimes", tokenList())()),
    defaultInterpolationScheme_("cell"),
    interpolationSchemes_(ITstream("interpolationSchemes", tokenList())()),
    particlesPerCell_(0),
    particleNumberControl_(),
    cloneAt_(),
    eliminateAt_(),
    kMin_("kMin", dimVelocity*dimVelocity, 100.0*SMALL),
    DNum_("DNum", dimless, 0.)
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mcSolution::read()
{
    if (regIOobject::read())
    {
        const dictionary& dict = solutionDict();

        if (dict.found("averagingTime"))
        {
            averagingTime_.value() = readScalar(dict.lookup("averagingTime"));
        }

        if (dict.found("relaxationTimes"))
        {
            relaxationTimes_ = dict.subDict("relaxationTimes");
        }

        if (relaxationTimes_.found("default"))
        {
            defaultRelaxationTime_.value() =
                readScalar(relaxationTimes_.lookup("default"));
        }

        if (dict.found("interpolationSchemes"))
        {
            interpolationSchemes_ = dict.subDict("interpolationSchemes");
        }

        if (interpolationSchemes_.found("default"))
        {
            interpolationSchemes_.lookup("default")
                >> defaultInterpolationScheme_;
        }

        particlesPerCell_ = readLabel(dict.lookup("particlesPerCell"));

        dict.lookup("particleNumberControl") >> particleNumberControl_;

        cloneAt_ = readScalar(dict.lookup("cloneAt"));
        if (cloneAt_ < 0 || !(cloneAt_ < 1.))
        {
            FatalErrorIn("mcSolution::read()")
                << "The value of " << dict.name() << "::cloneAt = " << cloneAt_
                << " must be in the range [0, 1)\n"
                << exit(FatalError);
        }

        eliminateAt_ = readScalar(dict.lookup("eliminateAt"));
        if (!(eliminateAt_ > 1))
        {
            FatalErrorIn("mcSolution::read()")
                << "The value of " << dict.name() << "::eliminateAt = "
                << eliminateAt_ << " must be > 1\n"
                << exit(FatalError);
        }

        if (dict.found("kMin"))
        {
            kMin_.value() = readScalar(dict.lookup("kMin"));
        }

        if (dict.found("DNum"))
        {
            DNum_.value() = readScalar(dict.lookup("DNum"));
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::mcSolution::readIfModified()
{
    if (regIOobject::readIfModified())
    {

        return read();
    }
    else
    {
        return false;
    }
}


const Foam::dictionary& Foam::mcSolution::solutionDict() const
{
    if (found("select"))
    {
        return subDict(word(lookup("select")));
    }
    else
    {
        return *this;
    }
}


Foam::dimensionedScalar
Foam::mcSolution::relaxationTime(const word& name) const
{
    if (relaxationTimes_.found(name))
    {
        return dimensionedScalar
            (
                name+"RelaxTime",
                dimTime,
                relaxationTimes_.lookup(name)
            );
    }
    return defaultRelaxationTime_;
}


Foam::word Foam::mcSolution::interpolationScheme(const word& name) const
{
    if (interpolationSchemes_.found(name))
    {
        return relaxationTimes_.lookup(name);
    }
    return defaultInterpolationScheme_;
}

// ************************************************************************* //
