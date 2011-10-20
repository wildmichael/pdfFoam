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
    Foam::mcParticleCloud

Description
    A Cloud of Monte Carlo Particles (mcParticle)

SourceFiles
    mcParticleCloudI.H
    mcParticleCloud.C
    mcParticleCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef mcParticleCloud_H
#define mcParticleCloud_H

#include "Cloud.H"
#include "mcParticle.H"
#include "IOdictionary.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class mcParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class mcParticleCloud
:
    public Cloud<mcParticle>
{
  // Private data
  
  const fvMesh& mesh_;
  
  const volVectorField& Ufv_;
  
  scalar dtCloud_;

  //- random number for Wiener process (diffusion) 
  Random random_;

  IOdictionary particleProperties_;

  //- statistical moments (mass, momentum, energy)
  //- needed to obtain UPdf, and TauPdf fields
  volScalarField M0;
  volVectorField M1;
  volSymmTensorField M2;

  // Private Member Functions

  //- Disallow default bitwise copy construct
  mcParticleCloud(const mcParticleCloud&);
  
  //- Disallow default bitwise assignment
  void operator=(const mcParticleCloud&);


public:

    // Constructors

        //- Construct given mesh
        mcParticleCloud
        (
            const fvMesh&,
            const volVectorField& U,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );


  // Member Functions

      // Access

      inline const fvMesh& mesh() const;
      inline scalar dtCloud();
      inline scalar e() const;
      inline scalar mu() const;

  
      Random &random() {return random_;}

        // Edit


  //- initial release of particles
  void initReleaseParticles();

  //- Update end position (for tracking)
  void updateParticleProperties();

  //- Evolve the particles
  void evolve();

  //- Report Info
  void info() const;
  void allParticlesInfo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mcParticleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
