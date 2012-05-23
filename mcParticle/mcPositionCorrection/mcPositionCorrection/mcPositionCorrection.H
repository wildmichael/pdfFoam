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
    Foam::mcPositionCorrection

Description
    Abstract base class for Monte Carlo position correction models

SourceFiles
    mcPositionCorrection.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcPositionCorrection_H
#define mcPositionCorrection_H

#include "mcModel.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class mcParticle;
class mcParticleCloud;

/*---------------------------------------------------------------------------*\
                     Class mcPositionCorrection Declaration
\*---------------------------------------------------------------------------*/

class mcPositionCorrection :
    public mcModel
{
    // Private Data

        //- Average of the inverse of the deltaCoeffs
        autoPtr<volScalarField> L_;

    // Private Member Functions

        void makeL();

        // Disallow default bitwise copy construct and assignment
        mcPositionCorrection(const mcPositionCorrection&);
        void operator=(const mcPositionCorrection&);

protected:

    // Protected Member Functions

        //- Average of the inverse of the deltaCoeffs
        const volScalarField& L();

public:

    //- Runtime type information
    TypeName("mcPositionCorrection");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mcPositionCorrection,
            mcPositionCorrection,
            (
                mcParticleCloud& cloud,
                const objectRegistry& db,
                const word& subDictName
            ),
            (cloud, db, subDictName)
        );

    // Constructors

        //- Construct from components
        mcPositionCorrection
        (
            mcParticleCloud& cloud,
            const objectRegistry& db,
            const word& subDictName
        );

    // Selectors

        //- Return a reference to the selected position correction
        static autoPtr<mcPositionCorrection> New
        (
            mcParticleCloud& cloud,
            const objectRegistry& db
        );

    // Member Functions

        using mcModel::correct;

        //- Apply the model to a single particle
        virtual void correct(mcParticle&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mcPositionCorrectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
