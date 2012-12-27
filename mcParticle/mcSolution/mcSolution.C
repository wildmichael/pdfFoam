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

#include "mcSolution.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcSolution::mcSolution(const objectRegistry& obr, const word& name)
:
    IOdictionary
    (
        IOobject
        (
            name,
            obr.time().system(),
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    CFL_(),
    averagingCoeff_(1000),
    defaultRelaxationTime_("defaultRelaxationTime", dimTime, GREAT),
    relaxationTimes_(ITstream("relaxationTimes", tokenList())()),
    minRelaxationFactor_(10),
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

        CFL_ = readScalar(dict.lookup("CFL"));

        if (dict.found("averagingCoeff"))
        {
            averagingCoeff_ = readScalar(dict.lookup("averagingCoeff"));
        }

        if (dict.found("relaxationTimes"))
        {
            relaxationTimes_ = dict.subDict("relaxationTimes");
        }

        bool bad = false;
        if (relaxationTimes_.found("default"))
        {
            token t(relaxationTimes_.lookup("default"));
            if (t.isWord())
            {
                word value = t.wordToken();
                if (value == "none")
                {
                    haveDefaultRelaxationTime_ = false;
                }
                else
                {
                    bad = true;
                }
            }
            else if (t.isNumber())
            {
                haveDefaultRelaxationTime_ = true;
                defaultRelaxationTime_.value() =
                    readScalar(relaxationTimes_.lookup("default"));
            }
            else
            {
                bad = true;
            }
            if (bad)
            {
                const fileName& entryName =
                    relaxationTimes_.lookupEntry
                    (
                        "default",
                        false,
                        false
                    ).name();
                FatalIOErrorIn("mcSolution::read()", relaxationTimes_)
                    << entryName
                    << " must either be a scalar or none\n"
                    << exit(FatalIOError);
            }
        }

        if (dict.found("minRelaxationFactor"))
        {
            minRelaxationFactor_ =
                readScalar(dict.lookup("minRelaxationFactor"));
        }

        if (dict.found("interpolationSchemes"))
        {
            interpolationSchemes_ = dict.subDict("interpolationSchemes");
        }

        if (interpolationSchemes_.found("default"))
        {
            interpolationSchemes_.lookup("default")
                >> defaultInterpolationScheme_;
            if (defaultInterpolationScheme_ == "none")
            {
                haveDefaultInterpolationScheme_ = false;
            }
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
    dimensionedScalar result = defaultRelaxationTime_;
    if (relaxationTimes_.found(name) || !haveDefaultRelaxationTime_)
    {
        result =
            dimensionedScalar
            (
                name+"RelaxTime",
                dimTime,
                relaxationTimes_.lookup(name)
            );
    }
    result.value() =
        max(result.value(), time().deltaT().value()*minRelaxationFactor_);
    return result;
}


Foam::word Foam::mcSolution::interpolationScheme(const word& name) const
{
    if (interpolationSchemes_.found(name) || !haveDefaultInterpolationScheme_)
    {
        return interpolationSchemes_.lookup(name);
    }
    return defaultInterpolationScheme_;
}

// ************************************************************************* //
