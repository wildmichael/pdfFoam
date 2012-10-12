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
    Foam::mcSteadyFlamelet

Description
    Steady flamelet reaction model

    Reads constant/flamelet/{z,chi,Rgas,T} (mixture fraction, scalar
    dissipation rate, specific gas constant and temperature) and additional
    fields specified in the "scalars" list. z must be uniformly spaced (up to
    0.1% of the first spacing) and monothonically increasing. chi must be
    monothonically increasing. If Nz is the number of entries in z and Nchi
    the number of entries in chi, then Rgas, T and all the fields specified
    in "scalars" must be lists of lists, the outer list of size Nchi, the inner
    of size Nz.

    The scalar dissipation rate chi is modelled using the model proposed by
    Poinsot and Veynante [1]

        chi = Cchi*Omega*favg(z''*z'')

    where Cchi is a model constant (here, 6.0 by default), Omega the particle
    turbulent frequency, rho the particle density, z'' the mixture fraction
    fluctuations and favg() the Favre-averaging operator.

    [1] T. Poinsot and D. Veynante, Theoretical and Numerical Combustion.
        Edwards, 2001.

SourceFiles
    mcSteadyFlamelet.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcSteadyFlamelet_H
#define mcSteadyFlamelet_H

#include "mcReactionModel.H"

#include "autoPtr.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class> class interpolation;

/*---------------------------------------------------------------------------*\
                       Class mcSteadyFlamelet Declaration
\*---------------------------------------------------------------------------*/

class mcSteadyFlamelet :
    public mcReactionModel
{

    // Private Data

        //- Names of the mixture-fraction and temperature fields
        const word zName_, TName_;

        //- Scalar dissipation-rate model constant
        const scalar Cchi_;

        //- Mixture-fraction and scalar dissipation-rate vectors
        const scalarIOList z_, chi_;

        //- Number of z and chi data points
        const label Nz_, Nchi_;

        //- Steps in z_
        scalar deltaZ_;

        //- Names of the additional interpolated fields
        wordList addedNames_;

        //- Interpolation tables
        PtrList<scalarListIOList> phi_;

        //- Indexes of the mixture fraction, scalar dissipation rate and
        // temperature properties
        const label zIdx_, chiIdx_, TIdx_;

        //- Indexes of the additional interpolated fields
        labelList addedIdx_;

        //- Number of native fields (temperature T and density rho for now)
        static const label nNativeFields_ = 2;

        //- The interpolator for the mixture fraction variance
        autoPtr<interpolation<scalar> > zVarInterp_;

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mcSteadyFlamelet(const mcSteadyFlamelet&);
        void operator=(const mcSteadyFlamelet&);

public:

    //- Runtime type information
    TypeName("mcSteadyFlamelet");

    // Constructors

        //- Construct from components
        mcSteadyFlamelet
        (
            mcParticleCloud& cloud,
            const objectRegistry& db,
            const word& subDictName
        );

    // Member Functions

        //- Update/prepare internal data
        virtual void updateInternals();

        //- Apply the reaction model to a single particle
        virtual void correct(mcParticle& p);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
