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

#include "mcThermo.H"

#include "IOdictionary.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(mcThermo, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mcThermo::calculate()
{
    mu_  = nu_ * rho_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcThermo::mcThermo(const fvMesh& mesh)
:
    basicThermo(mesh),

    mesh_(mesh),

    cloudP_(0),

    rho_
    (
        IOobject
        (
            lookupOrDefault<word>("rhoName", "rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    nu_
    (
        IOdictionary
        (
            IOobject
            (
                "transportProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).lookup("nu")
    ),

    nFVSubCycles_(lookupOrDefault<label>("nFVSubCycles", 1)),
    nPDFSubCycles_(lookupOrDefault<label>("nPDFSubCycles", 1))
{
    calculate();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mcThermo::~mcThermo()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcThermo::correct()
{
    if (debug)
    {
        Info<< "entering mcThermo::correct()" << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "exiting mcThermo::correct()" << endl;
    }
}


void Foam::mcThermo::evolve()
{
    Time& runTime = const_cast<Time&>(mesh_.time());
    if ((runTime.timeIndex()-runTime.startTimeIndex()) % (nPDFSubCycles_ + nFVSubCycles_) + 1 > nFVSubCycles_)
    {
        if (debug)
        {
            Info<< "executing mcThermo::cloudP_().evolve()" << endl;
            Info<< tab << "cloudP_().size() = ";
            if (cloudP_.valid())
               Info << cloudP_().size();
            else
               Info << "empty";
           Info << endl;
        }

        // Instantiate the cloud on first call
        if (!cloudP_.valid())
        {
            cloudP_.reset
            (
                new mcParticleCloud
                (
                    mesh_,
                    subDict("cloudProperties"),
                    lookupOrDefault<word>("cloudName", "mcThermoCloud"),
                    0, 0, &rho_
                )
            );
        }
        for (label i=0; i < nPDFSubCycles_; i++)
        {
            runTime++;
            cloudP_().evolve();
        }

        if (debug)
        {
            Info<< "done executing mcThermo::cloudP_().evolve()" << endl;
            Info<< tab << "cloudP_().size() = " << cloudP_().size() << endl;
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::mcThermo::rho() const
{
    return rho_;
}


Foam::volScalarField& Foam::mcThermo::rho()
{
    return rho_;
}

// ************************************************************************* //
