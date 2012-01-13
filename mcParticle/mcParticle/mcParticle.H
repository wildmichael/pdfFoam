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
#include "interpolationCellPointFace.H"
#include "gradInterpolationConstantTet.H"
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

        //- particle correction velocity
        vector Ucorrection_;

        //- particle tracking velocity leading to destination position
        vector Utracking_;

       //- fluid velocity @ particle location, updated at the end of each Lagrangian dt
       // (same time with Updf & UParticle)
       vector UFap_;

       //- turbulent frequency of particle
       scalar Omega_;

       //- density
       scalar rho_;

       //- separation (only for ghost cells)
       vector shift_;

       //- particle Courant number
       scalar Co_;

       //- Velocity of the "moving boundary" the particle bounced of
       vector reflectionBoundaryVelocity_;

       // ==== PUT ALL label AND OTHER 4-BYTE TYPES AFTER THIS === //

       //- Whether this particle is a ghost
       label ghost_;

       //- Number of tracking steps during single time-step
       label nSteps_;

       // ==== PUT ALL char AND OTHER 2-BYTE TYPES AFTER THIS === //

       // ==== PUT ALL bool AND OTHER 1-BYTE TYPES AFTER THIS === //

       //- Whether this particle is located on an inlet boundary
       // and must be assigned a random stepFraction
       bool isOnInletBoundary_;

       //- Whether this particle has been reflected at an mcOpenBoundary
       bool reflectedAtOpenBoundary_;

       // ==== PUT ALL NON-POD DATA AFTER THIS === //

       //- scalar properties
       scalarField Phi_;

public:

    friend class Cloud<mcParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackData:
        public Particle<mcParticle>::trackData
    {

        // Interpolators for continuous phase fields

            const interpolationCellPointFace<scalar>&   rhoInterp_;
            const interpolationCellPointFace<vector>&   UInterp_;
            const gradInterpolationConstantTet<scalar>& gradPInterp_;
            const interpolationCellPointFace<scalar>&   kInterp_;
            const interpolationCellPointFace<vector>&   diffUInterp_;


    public:

        // Constructors

        trackData
        (
            mcParticleCloud& mcpc,
            const interpolationCellPointFace<scalar>&   rhoInterp,
            const interpolationCellPointFace<vector>&   UInterp,
            const gradInterpolationConstantTet<scalar>& gradPInterp,
            const interpolationCellPointFace<scalar>&   kInterp,
            const interpolationCellPointFace<vector>&   diffUInterp
        );

        // Member functions

        inline const interpolationCellPointFace<scalar>& rhoInterp() const;

        inline const interpolationCellPointFace<vector>& UInterp() const;

        inline const gradInterpolationConstantTet<scalar>& gradPInterp() const;

        inline const interpolationCellPointFace<scalar>& kInterp() const;

        inline const interpolationCellPointFace<vector>& diffUInterp() const;
    };


    // Constructors

        //- Construct from components
        mcParticle
        (
            const mcParticleCloud& c,
            const vector& position,
            const label celli,
            const scalar m,
            const vector& Updf,
            const vector& UParticle,
            const vector& UFap,
            const scalarField& Phi,
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

            //- Return total particle velocity;
            inline vector& UParticle();

            //- Return particle correction velocity
            inline const vector& Ucorrection() const;

            //- Return particle correction velocity
            inline vector& Ucorrection();

            //- Return velocity with which the particle is being tracked
            inline const vector& Utracking() const;

            //- Fluid velocity @ particle location
            inline const vector& UFap() const;

            //- Turbulent frequency
            inline scalar Omega() const;

            //- Turbulent frequency
            inline scalar& Omega();

            //- scalar properties
            inline const scalarField& Phi() const;

            //- scalar properties
            inline scalarField& Phi();

            //- density
            inline scalar rho() const;

            //- density
            inline scalar& rho();

            //- whether this particle is ghost
            inline label ghost() const;

            //- whether this particle is ghost
            inline label& ghost();

            //- number of tracking steps during single time step
            inline label nSteps() const;

            //- number of tracking steps during single time step
            inline label& nSteps();

            //- shift for ghost particles
            inline const vector& shift() const;

            //- shift for ghost particles
            inline vector& shift();

            //- particle Courant number
            inline scalar Co() const;

            //- Velocity of the "moving boundary" the particle bounced of
            inline const vector& reflectionBoundaryVelocity() const;

            //- Velocity of the "moving boundary" the particle bounced of
            inline vector& reflectionBoundaryVelocity() ;

            //- whether this particle is located on a inlet boundary
            inline bool isOnInletBoundary() const;

            //- whether this particle is located on a inlet boundary
            inline bool& isOnInletBoundary();

            //- Whether this particle has been reflected at an mcOpenBoundary
            inline bool reflectedAtOpenBoundary() const;

            //- Whether this particle has been reflected at an mcOpenBoundary
            inline bool& reflectedAtOpenBoundary();


        // Tracking

            //- Move
            bool move(trackData&);


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch& patch,
                trackData&       td,
                const label      patchI
            );

            //- Overridable function to handle the particle hitting a patch
            //  without trackData.
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch& patch,
                int&             td,
                const label      patchI
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                mcParticle::trackData& td
            )
            {}

            //- Overridable function to handle the particle hitting a polyPatch
            //- without trackData
            void hitPatch
            (
                const polyPatch&,
                int&
            )
            {}


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

            //- Info string about this particle and its properties
            string info() const;


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
