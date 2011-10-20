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

SourceFiles
    mcRASOmegaModel.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcRASOmegaModel_H
#define mcRASOmegaModel_H

#include "mcOmegaModel.H"

#include "autoPtr.H"
#include "interpolationCellPoint.H"

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

        //- Interpolator for Omega
        autoPtr<interpolationCellPoint<scalar> > OmegaInterp_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcRASOmegaModel(const mcRASOmegaModel&);
        void operator=(const mcRASOmegaModel&);

        //- Setup the interpolator
        void setupInterpolator(const mcParticleCloud& cloud);

public:

    //- Runtime type information
    TypeName("mcRASOmegaModel");

    // Constructors

        //- Construct from components
        mcRASOmegaModel
        (
            const dictionary& mcRASOmegaModelDict
        );


    // Member Functions

        //- Apply the Omega model to the particle cloud
        virtual void correct(mcParticleCloud& cloud);

        //- Apply the Omega model to a single particle
        virtual void correct
        (
            mcParticleCloud& cloud,
            mcParticle& p,
            bool prepare=true
        );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
