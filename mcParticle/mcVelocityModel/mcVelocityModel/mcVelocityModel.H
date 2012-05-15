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

Class
    Foam::mcVelocityModel

Description
    Abstract base class for Monte Carlo particle velocity models

SourceFiles
    mcVelocityModel.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcVelocityModel_H
#define mcVelocityModel_H

#include "mcModel.H"

#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "mcParticle.H"
#include "runTimeSelectionTables.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mcVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class mcVelocityModel :
    public mcModel
{
    // Private Member Data

        //- The cloud
        mcParticleCloud& cloud_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcVelocityModel(const mcVelocityModel&);
        void operator=(const mcVelocityModel&);

protected:

    // Protected Member Functions

        //- Update the internal fields
        virtual void setupInternals();

        //- Return the cloud
        mcParticleCloud& cloud() {return cloud_;}

        //- Return the cloud
        const mcParticleCloud& cloud() const {return cloud_;}

public:

    //- Runtime type information
    TypeName("mcVelocityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mcVelocityModel,
            mcVelocityModel,
            (
                const Foam::objectRegistry& db,
                const dictionary& parentDict,
                const dictionary& dict,
                mcParticleCloud& cloud
            ),
            (db, parentDict, dict, cloud)
        );

    // Constructors

        //- Construct from components
        mcVelocityModel
        (
            const Foam::objectRegistry& db,
            const dictionary& parentDict,
            const dictionary& mcVelocityModelDict,
            mcParticleCloud& cloud
        );


    // Selectors

        //- Return a reference to the selected velocity model
        static autoPtr<mcVelocityModel> New
        (
            const Foam::objectRegistry& db,
            const dictionary& dict,
            mcParticleCloud& cloud
        );

    //- Destructor
        virtual ~mcVelocityModel()
        {}

    // Member Functions

        //- Apply the velocity model to the particle cloud
        virtual void correct();

        //- Apply the velocity model to a single particle
        virtual void correct
        (
            mcParticle& p,
            bool prepare=true
        ) = 0;

        //- Apply the velocity model to the particle cloud
        virtual void correct(mcParticleCloud&) {correct();}

        //- The mean turbulence frequency
        virtual tmp<volScalarField> omega() const = 0;

        //- The diffusive coefficient based on the velocity model
        // Rationale:
        //     dX = U*dt
        //     dU = -1/tau*U*dt + sqrt(2*D)*dW
        //   for tau -> 0
        //     dX = sqrt(2*D*sqr(tau))*dW
        //   where
        //     tau = (0.5 + 0.75*C0)*omega
        //     D = k*omega*C0/2
        //   it follows
        //     dX = (sqrt(omega*mut/rho*Cdiff)
        //   where Cdiff is the diffusive coefficient
        //     Cdiff = C0/Cmu*sqr(0.5 + 0.75*C0)
        virtual dimensionedScalar Cdiff() const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
