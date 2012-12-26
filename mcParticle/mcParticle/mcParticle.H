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

#if FOAM_HEX_VERSION < 0x200
#include "Particle.H"
#else
#include "particle.H"
#endif
#include "IOstream.H"
#include "autoPtr.H"
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
#if FOAM_HEX_VERSION < 0x200
    public Particle<mcParticle>
#else
    public particle
#endif
{
    // Private types

#if FOAM_HEX_VERSION < 0x200
        typedef Particle<mcParticle> base;
#else
        typedef particle base;
#endif

        struct BoundaryOfDataMembers {} beginOfDataMembers_;

        // ==== NO DATA MEMBERS BEFORE THIS === //

    // Private member data

        //- old particle position
        point positionOld_;

        //- old cell
        label celliOld_;

        //- old face
        label faceiOld_;

        //- old processor
        label procOld_;

        //- mass (statistical weight)
        scalar m_;

        //- particle velocity (mean + flutuation) (primary variable)
        vector UParticle_;

        //- old particle velocity
        vector UParticleOld_;

        //- particle correction velocity
        vector Ucorrection_;

        //- particle tracking velocity leading to destination position
        vector Utracking_;

        //- turbulent frequency of particle
        scalar Omega_;

        //- density
        scalar rho_;

        //- eta
        scalar eta_;

        //- separation (only for ghost cells)
        vector shift_;

        //- particle Courant number
        scalar Co_;

        //- Velocity of the "moving boundary" the particle bounced of
        vector reflectionBoundaryVelocity_;

        //- Whether this particle is a ghost
        label ghost_;

        //- Number of tracking steps during single time-step
        label nSteps_;


        //- Whether this particle is located on an inlet boundary
        // and must be assigned a random stepFraction
        bool isOnInletBoundary_;

        //- Whether this particle has been reflected at any type of patch
        bool reflected_;

        //- Whether this particle has been reflected at an mcOpenBoundary
        bool reflectedAtOpenBoundary_;

        // ==== NO POD DATA MEMBERS AFTER THIS === //

        BoundaryOfDataMembers endOfDataMembers_;

        // ==== PUT ALL NON-POD DATA AFTER THIS === //

        //- scalar properties
        scalarField Phi_;

public:

#if FOAM_HEX_VERSION >= 0x200
    //- Runtime type information
    TypeName("mcParticle");
#endif

    friend class Cloud<mcParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackData:
#if FOAM_HEX_VERSION < 0x200
        public Particle<mcParticle>::trackData
#else
        public particle::TrackingData<mcParticleCloud>
#endif
    {
#if FOAM_HEX_VERSION < 0x200
        typedef Particle<mcParticle>::trackData base;
#else
        typedef particle::TrackingData<mcParticleCloud> base;
#endif
        //- The cloud
        mcParticleCloud& cloud_;

#if FOAM_HEX_VERSION < 0x200
        //- The time span over which to track the particle (deltaT)
        scalar trackTime_;

        //- The global deltaT
        scalar deltaT_;
#endif

    public:

        // Constructors

        trackData(mcParticleCloud& mcpc, scalar trackTime);

        // Access

            //- Return a reference to the cloud
            mcParticleCloud& cloud() { return cloud_; }

#if FOAM_HEX_VERSION < 0x200

            //- The time span over which to track the particle (deltaT)
            scalar trackTime() const { return trackTime_; }

            //- The global time step
            scalar deltaT() const { return deltaT_; }
#endif
    };


    // Constructors

        //- Construct from components
        mcParticle
        (
            const mcParticleCloud& c,
            const vector& position,
            const label celli,
            const scalar m,
            const vector& UParticle,
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

#if FOAM_HEX_VERSION < 0x200
        //- Construct and return a clone
        autoPtr<mcParticle> clone() const
        {
            return autoPtr<mcParticle>(new mcParticle(*this));
        }
#else
        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new mcParticle(*this));
        }


        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<mcParticle> operator()(Istream& is) const
            {
                return autoPtr<mcParticle>
                (
                    new mcParticle(mesh_, is, true)
                );
            }
        };
#endif


    // Member Functions

        // Access

            //- Face index if the particle is on a face otherwise -1
            using base::face;

#if FOAM_HEX_VERSION < 0x200
            //- Face index if the particle is on a face otherwise -1
            inline label& face();
#endif

            //- Return the old position
            inline const point& positionOld() const;

            //- Return the old position
            inline point& positionOld();

            //- Return the old cell
            inline label cellOld() const;

            //- Return the old cell
            inline label& cellOld();

            //- Return the old face
            inline label faceOld() const;

            //- Return the old face
            inline label& faceOld();

            //- Return the old processor
            inline label procOld() const;

            //- Return the old processor
            inline label& procOld();

            //- Return mass
            inline scalar m() const;

            inline scalar& m();

           // inline const symmTensor& Taupdf() const;
           // inline symmTensor& Taupdf();

            //- Return total particle velocity;
            inline const vector& UParticle() const;

            //- Return total particle velocity;
            inline vector& UParticle();

            //- Return old particle velocity
            inline const vector& UParticleOld() const;

            //- Return old particle velocity
            inline vector& UParticleOld();

            //- Return particle correction velocity
            inline const vector& Ucorrection() const;

            //- Return particle correction velocity
            inline vector& Ucorrection();

            //- Return velocity with which the particle is being tracked
            inline const vector& Utracking() const;

            //- Return velocity with which the particle is being tracked
            inline vector& Utracking();

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

            //- time stepping adjusting parameter
            inline scalar eta() const;

            //- time stepping adjusting parameter
            inline scalar& eta();

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

            //- particle Courant number
            inline scalar& Co();

            //- Velocity of the "moving boundary" the particle bounced of
            inline const vector& reflectionBoundaryVelocity() const;

            //- Velocity of the "moving boundary" the particle bounced of
            inline vector& reflectionBoundaryVelocity() ;

            //- whether this particle is located on a inlet boundary
            inline bool isOnInletBoundary() const;

            //- whether this particle is located on a inlet boundary
            inline bool& isOnInletBoundary();

            //- Whether this particle has been reflected at any type of patch
            inline bool reflected() const;

            //- Whether this particle has been reflected at any type of patch
            inline bool& reflected();

            //- Whether this particle has been reflected at an mcOpenBoundary
            inline bool reflectedAtOpenBoundary() const;

            //- Whether this particle has been reflected at an mcOpenBoundary
            inline bool& reflectedAtOpenBoundary();


        // Tracking

            //- Move
#if FOAM_HEX_VERSION < 0x200
            bool move(trackData&);
#else
            bool move(trackData& td, const scalar trackTime);
#endif


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch&  patch,
                trackData&        td,
#if FOAM_HEX_VERSION < 0x200
                const label       patchI
#else
                const label       patchI,
                const scalar      trackFraction,
                const tetIndices& tetIs
#endif
            );

#if FOAM_HEX_VERSION < 0x200
            //- Overridable function to handle the particle hitting a patch
            //  without trackData.
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch& patch,
                int&             td,
                const label      patchI
            );
#endif

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                trackData& td
            )
            {}

#if FOAM_HEX_VERSION < 0x200
            //- Overridable function to handle the particle hitting a polyPatch
            //- without trackData
            void hitPatch
            (
                const polyPatch&,
                int&
            )
            {}
#endif


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
