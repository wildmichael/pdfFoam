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
    Foam::basicPdfThermo

Description
    Density from ideal gas law: rho=p/(R*T)

SourceFiles
    basicPdfThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicPdfThermo_H
#define basicPdfThermo_H

#include "basicRhoThermo.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class basicPdfThermo Declaration
\*---------------------------------------------------------------------------*/

class basicPdfThermo
:
    public basicRhoThermo
{
    // Private data

        //- Particle cloud
        mcParticleCloud cloud_;

        //- Viscosity
        dimensionedScalar nu_;

        //- Evolution frequency
        label evolutionFrequency_;

    // Private member functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        basicPdfThermo(const basicPdfThermo&);


public:

    //- Runtime type information
    TypeName("basicPdfThermo");


    // Constructors

        //- Construct from mesh
        basicPdfThermo(const fvMesh& mesh);


    //- Destructor
    virtual ~basicPdfThermo();


    // Member functions

        //- Evolve the cloud
        void evolve();

        //- Update properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
