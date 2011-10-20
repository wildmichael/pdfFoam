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

        //- Reference to mean velocity
        vector Updf_;
  
        //- Fluctuation Velocity
        vector u_;

       //- current destination
        vector endPosition_;

       //- current time step
        scalar dt_;

public:

    friend class Cloud<mcParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackData
    {
        //- Reference to the cloud containing this particle
        mcParticleCloud& mcpc_;

        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& rhoInterp_;
            const interpolationCellPoint<vector>& UInterp_;
            const interpolationCellPoint<scalar>& kInterp_;


    public:

      bool switchProcessor;
      bool keepParticle;
      
      // Constructors
      
      inline trackData
      (
       mcParticleCloud& mcpc,
       const interpolationCellPoint<scalar>& rhoInterp,
       const interpolationCellPoint<vector>& UInterp,
       const interpolationCellPoint<scalar>& kInterp
       );

      // Member functions

      inline mcParticleCloud& mcpc();
      
      inline const interpolationCellPoint<scalar>& rhoInterp() const;
      
      inline const interpolationCellPoint<vector>& UInterp() const;
      
      inline const interpolationCellPoint<scalar>& kInterp() const;
      
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
            const vector& u,
            const scalar dt
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

            //- Return velocity PDF
            inline const vector& Updf() const;

            //- Return fluctuation
            inline const vector& u() const;

            //- Return total particle velocity;
            inline const vector UParticle();
            
            //- particle time step
            inline scalar dt();

           //- Return current destination
            inline vector& endPosition();

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
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mcParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
