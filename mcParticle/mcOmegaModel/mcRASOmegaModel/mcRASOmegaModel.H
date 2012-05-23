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
    Foam::mcRASOmegaModel

Description
    Omega model based on data from the RAS model

    The coefficient
    system/mcSolution::kMin can be used to put a lower bound on the TKE when
    computing the turbulence frequency from epsilon/k.

SourceFiles
    mcRASOmegaModel.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcRASOmegaModel_H
#define mcRASOmegaModel_H

#include "mcOmegaModel.H"

#include "autoPtr.H"
#include "interpolationCellPointFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;

/*---------------------------------------------------------------------------*\
                       Class mcRASOmegaModel Declaration
\*---------------------------------------------------------------------------*/

class mcRASOmegaModel :
    public mcOmegaModel
{
    // Private Data

        //- Omega field
        autoPtr<volScalarField> Omega_;

        //- Interpolator for Omega
        autoPtr<interpolationCellPointFace<scalar> > OmegaInterp_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcRASOmegaModel(const mcRASOmegaModel&);
        void operator=(const mcRASOmegaModel&);

public:

    //- Runtime type information
    TypeName("mcRASOmegaModel");

    // Constructors

        //- Construct from components
        mcRASOmegaModel
        (
            mcParticleCloud& cloud,
            const objectRegistry& db,
            const word& subDictName
        );


    // Member Functions

        //- Setup the interpolator
        void updateInternals();

        //- Apply the Omega model to a single particle
        virtual void correct(mcParticle& p);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
