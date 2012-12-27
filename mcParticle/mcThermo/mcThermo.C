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
    )
{
    calculate();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mcThermo::~mcThermo()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcThermo::createCloud()
{
    // Instantiate the cloud on first call
    if (!cloudP_.valid())
    {
        word cloudName = lookupOrDefault<word>("cloudName", "mcThermoCloud");
        if (debug)
        {
            Info<< "\nCreating mcParticleCloud " << cloudName << nl << endl;
        }
        cloudP_.reset
        (
            new mcParticleCloud
            (
                mesh_,
                subDict(cloudName+"Properties"),
                cloudName,
                0, 0, 0, &rho_
            )
        );
    }
    else
    {
        FatalErrorIn("mcThermo::createCloud()")
            << "Cloud already created!\n"
            << abort(FatalError);
    }
}


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


Foam::scalar Foam::mcThermo::evolve()
{
    Info<< "Evolving Monte Carlo particle cloud " << cloudP_().name() << nl
        << endl;
    return cloudP_().evolve();
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
