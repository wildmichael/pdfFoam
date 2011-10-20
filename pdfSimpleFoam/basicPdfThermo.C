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

#include "basicPdfThermo.H"

#include "IOdictionary.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basicPdfThermo, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::basicPdfThermo::calculate()
{
    mu_  = nu_ * rho_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicPdfThermo::basicPdfThermo(const fvMesh& mesh)
:
    basicRhoThermo(mesh),

    cloud_
    (
        mesh,
        subDict("cloudProperties"),
        0,
        0,
        0,
        lookupOrDefault<word>("cloudName", "pdfThermoCloud")
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

    evolutionFrequency_(lookupOrDefault<label>("evolutionFrequency", 1))
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicPdfThermo::~basicPdfThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::basicPdfThermo::correct()
{
    if (debug)
    {
        Info<< "entering basicPdfThermo::correct()" << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "exiting basicPdfThermo::correct()" << endl;
    }
}

void Foam::basicPdfThermo::evolve()
{
    if (!(rho_.mesh().time().timeIndex() % evolutionFrequency_))
    {
        if (debug)
        {
            Info<< "executing basicPdfThermod::cloud_.evolve()" << endl;
            Info<< tab << "cloud_.size() = " << cloud_.size() << endl;
        }

        cloud_.evolve();

        if (debug)
        {
            Info<< "done executing basicPdfThermod::cloud_.evolve()" << endl;
            Info<< tab << "cloud_.size() = " << cloud_.size() << endl;
        }
    }
}


// ************************************************************************* //
