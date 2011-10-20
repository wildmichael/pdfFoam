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
    Foam::mcParticle

Description
    Monte Carlo Particle used in PDF method.

SourceFiles
    mcParticleI.H
    mcParticle.C
    mcParticleIO.C

Author: Heng Xiao with Patrick Jenny, 2010-08-29

\*---------------------------------------------------------------------------*/

#ifndef mcParticle_H
#define mcParticle_H

#include "Particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"
#include "meshTools.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mcParticleCloud;

/*---------------------------------------------------------------------------*\
                           Class mcParticle Declaration
\*---------------------------------------------------------------------------*/

class mcParticle
:
    public Particle<mcParticle>
{

    // Private member data

        //- mass (statistical weight)
        scalar m_;

        //- mean velocity based on particles
        vector Updf_;

        //- particle velocity (mean + flutuation) (primary variable)
        vector UParticle_;

       //- fluid velocity @ particle location, updated at the end of each Lagrangian dt
       // (same time with Updf & UParticle)
       vector UFap_;

       //- density
       scalar rho_;

       //- current time step
       scalar dt_;

       //- separation (only for ghost cells)
       vector shift_;

       // ==== PUT ALL label AND OTHER 4-BYTE TYPES AFTER THIS === //

       //- Whether this particle is a ghost
       label ghost_;

       // ==== PUT ALL NON-POD DATA AFTER THIS === //

       //- scalar properties
       scalarField Phi_;

protected:

       //- update the thermodynamic data
       void updateThermo();

public:

    friend class Cloud<mcParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackData:
        public Particle<mcParticle>::trackData
    {

        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& rhoInterp_;
            const interpolationCellPoint<vector>& UInterp_;
            const interpolationCellPoint<vector>& gradPInterp_;
            const interpolationCellPoint<scalar>& kInterp_;
            const interpolationCellPoint<scalar>& epsilonInterp_;
            const PtrList<interpolationCellPoint<scalar> >& PhiInterp_;
            const interpolationCellPoint<vector>& gradRhoInterp_;
            const interpolationCellPoint<vector>& diffUInterp_;


    public:


      // Constructors

      trackData
      (
       mcParticleCloud& mcpc,
       const interpolationCellPoint<scalar>& rhoInterp,
       const interpolationCellPoint<vector>& UInterp,
       const interpolationCellPoint<vector>& gradPInterp,
       const interpolationCellPoint<scalar>& kInterp,
       const interpolationCellPoint<scalar>& epsilonInterp,
       const PtrList<interpolationCellPoint<scalar> >& PhiInterp,
       const interpolationCellPoint<vector>& gradRhoInterp,
       const interpolationCellPoint<vector>& diffUInterp
       );

      // Member functions


      inline const interpolationCellPoint<scalar>& rhoInterp() const;

      inline const interpolationCellPoint<vector>& UInterp() const;

      inline const interpolationCellPoint<vector>& gradPInterp() const;

      inline const interpolationCellPoint<scalar>& kInterp() const;

      inline const interpolationCellPoint<scalar>& epsilonInterp() const;

      inline const PtrList<interpolationCellPoint<scalar> >& PhiInterp() const;

      inline const interpolationCellPoint<vector>& gradRhoInterp() const;

      inline const interpolationCellPoint<vector>& diffUInterp() const;
    };


    // Constructors

        //- Construct from components
        inline mcParticle
        (
            const Cloud<mcParticle>& c,
            const vector& position,
            const label celli,
            const scalar m,
            const vector& Updf,
            const vector& UParticle,
            const vector& UFap,
            const scalarField& Phi,
            const scalar dt,
            const vector& shift = vector::zero,
            const label  ghost = 0

        );

        //- Construct from Istream
        mcParticle
        (
            const Cloud<mcParticle>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<mcParticle> clone() const
        {
            return autoPtr<mcParticle>(new mcParticle(*this));
        }


    // Member Functions

        // Access

            //- Return mass
            inline scalar m() const;

            inline scalar& m();

            //- Return pdf velocity (particle velocity in this cell)
            inline const vector& Updf() const;

            inline vector& Updf();

           // inline const symmTensor& Taupdf() const;
           // inline symmTensor& Taupdf();

            //- Return total particle velocity;
            inline const vector& UParticle() const;

            //- Fluid velocity @ particle location
            inline const vector& UFap() const;

            //- scalar properties
            inline const scalarField& Phi() const;

            //- scalar properties
            inline scalarField& Phi();

            //- density
            inline scalar rho() const;

            //- density
            inline scalar& rho();

            //- particle time step
            inline scalar& dt();

            //- whether this particle is ghost
            inline label ghost() const;

            //- whether this particle is ghost
            inline label& ghost();

            //- shift for ghost particles
            inline const vector& shift() const;

            //- shift for ghost particles
            inline vector& shift();


        // Tracking

            //- Move
            bool move(trackData&);


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch&,
                mcParticle::trackData& td,
                const label patchI
            );

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions without trackData
            bool hitPatch
            (
                const polyPatch& p,
                int& td,
                const label patchI
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                mcParticle::trackData& td
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch without trackData
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                const wallPolyPatch&,
                mcParticle::trackData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            //- without trackData
            void hitWallPatch
            (
                const wallPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                mcParticle::trackData& td
            );

            //- Overridable function to handle the particle hitting a polyPatch
            //- without trackData
            void hitPatch
            (
                const polyPatch&,
                int&
            );

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            void transformProperties
            (
                const tensor& T
            );

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            void transformProperties
            (
                const vector& separation
            );


    // I-O

        static void readFields(Cloud<mcParticle>& c);

        static void writeFields(const Cloud<mcParticle>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const mcParticle&);
};


template<>
inline bool contiguous<mcParticle>()
{
    // contains scalar field
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mcParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
