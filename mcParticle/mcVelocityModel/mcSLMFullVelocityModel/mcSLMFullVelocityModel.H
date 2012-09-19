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
    Foam::mcSLMFullVelocityModel

Description
    Velocity model based on the simplified Langevin equation for the full
    velocity vector.

    See P. Jenny et al., "A Hybrid Algorithm for the Joint PDF Equation of
    Turbulent Reactive Flows", JCP 166 218-252,
    http://dx.doi.org/10.1006/jcph.2000.6646.

    This model adds correction drift terms for the mean velocity
    (@c system/mcSolution::relaxationTimes::U) and the turbulent kinetic energy
    (@c system/mcSolution::relaxationTimes::k). The SLM constant is controlled
    by the @c C0 dictionary entry (defaults to 2.1). Do not set @c C1; it is
    here for testing purposes only, switching off dissipation if set to 0.

SourceFiles
    mcSLMFullVelocityModel.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcSLMFullVelocityModel_H
#define mcSLMFullVelocityModel_H

#include "mcVelocityModel.H"

#include "autoPtr.H"
#include "gradInterpolationConstantTet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;
template<class> class interpolation;

/*---------------------------------------------------------------------------*\
                    Class mcSLMFullVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class mcSLMFullVelocityModel :
    public mcVelocityModel
{
    // Private Data

        //- The FV pressure
        const volScalarField& pfv_;

        //- The physical pressure
        volScalarField p_;

        //- The velocity difference
        volVectorField diffU_;

        //- The TKE difference
        scalarField diffk_;

        //- The decomposed mesh
        tetFacePointCellDecomposition<richTetPointRef> tetDecomp_;

        //- Interpolators
        autoPtr<gradInterpolationConstantTet<scalar> > gradPInterp_;
        autoPtr<interpolation<scalar> >   kcPdfInterp_;
        autoPtr<interpolation<scalar> >   kInterp_;
        autoPtr<interpolation<scalar> >   rhoInterp_;
        autoPtr<interpolation<vector> >   diffUInterp_;
        autoPtr<interpolation<vector> >   UInterp_;

        // Model parameters

            //- Simplified Langevin Model constant
            scalar C0_;

            //- Coefficient for the dissipation in the SLM
            scalar C1_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcSLMFullVelocityModel(const mcSLMFullVelocityModel&);
        void operator=(const mcSLMFullVelocityModel&);

public:

    //- Runtime type information
    TypeName("mcSLMFullVelocityModel");

    // Constructors

        //- Construct from components
        mcSLMFullVelocityModel
        (
            mcParticleCloud& cloud,
            const objectRegistry& db,
            const word& subDictName
        );


    // Member Functions

        //- Update/prepare internal data
        virtual void updateInternals();

        //- Apply the Velocity model to a single particle
        virtual void correct(mcParticle& p);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
