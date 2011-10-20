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
#include "labelIOField.H"
#include "SortableList.H"
#include "vectorList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

  typedef SortableList<mcParticle*> mcParticleList;
  typedef List <mcParticleList> mcParticleListList;

/*---------------------------------------------------------------------------*\
                           Class mcParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class mcParticleCloud
:
    public Cloud<mcParticle>
{
  
  enum cellPopStatus {EMPTY, NORMAL, TOOFEW, TOOMANY};
  //- scalar comparison function used for sorting
       class less
       {
       public:
    
         inline bool operator()
         (
           mcParticle* one,
           mcParticle* two
          ) const;
       };

      class more
      {
      public:
        
        inline bool operator()
        (
         mcParticle* one,
         mcParticle* two
         ) const;
      };
  

  // Private data
  bool  debug_;
  const fvMesh& mesh_;
  
  const volVectorField& Ufv_;
  const volScalarField& rhofv_;
  const volScalarField& kfv_;
  const volScalarField& psifv_;
  IOdictionary particleProperties_;
  
  scalar dtCloud_;
  scalar AvgTimeScale_;


  //- random number for Wiener process (diffusion) 
  Random random_;
  //- Number of particles per cell (initially)
  label Npc_;
  //- Cell size
  const label Nc_;
  // How many particle existed in history (including living ones)
  // only include particles generated in this run.
  scalar histNp_; 

  //- statistical moments (mass, momentum, energy)
  //- needed to obtain Updf (and later Taupdf) fields
  dimensionedScalar SMALL_MASS;

  labelIOField   PaNIC_; // Particle # in cell
  mcParticleListList  cellParticleAddr_;
  volScalarField M0_;
  volVectorField M1_;
  volScalarField Mpsi1_;
  volSymmTensorField M2_;
  scalarField instantM0_;

  labelListList   ghostCellLayers_;
  labelList       ghostPatchId_;
  labelListList   ghostFaceLayers_;
  HashSet<label>  ghostCellHash_;
  HashSet<label>  ghostFaceHash_;

  // Mesh based velocity and Reynolds stress, obtained from particle
  // quantities.
  volScalarField rhocPdf_;
  volVectorField UcPdf_;
  volScalarField psicPdf_;
  volSymmTensorField TaucPdf_;

  //- constants

  dimensionedScalar coeffRhoCorr_;
  dimensionedScalar coeffUCorr_;

  // Private Member Functions

  //- ensure moments M0, M1, M2 are cnosistently read
  void checkMoments();

  //- update moments and the quantities remembered by particles
  void updateCloudPDF(scalar existWt);
  void updateParticlePDF();

  //Particle number control actions
  void particleNumberControl();
  void cloneParticles(label celli);
  void clusterParticles(label celli);

  // Enforce in/out flow BCs by populating ghost cells
  void findGhostLayers();
  void populateGhostCells();

  // Boundary Fields for cloudPDF
  void setCloudPDFBoundary();
  void correctCloudPDFBoundaryConditions();

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
            const volScalarField& rho,
            const volScalarField& k,
            const volScalarField& psi,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );


  // Member Functions

      // Access

      inline const fvMesh& mesh() const;
      inline scalar e() const;
      inline scalar mu() const;

      inline scalar dtCloud();
      inline label Npc();
  
      Random &random() {return random_;}

        // Edit


  //- initial release of particles
  void initReleaseParticles();

  // Randomly generate N particles in celli, with provided cell-based
  // values and the scale of velocity fluctuation
  void particleGenInCell
  (
   label celli, 
   label N, 
   scalarList masses, 
   vector Updf, 
   vectorList uscales,
   scalar psi
   );

  // with the same properties
  void particleGenInCell
  (
   label celli, 
   label N, 
   scalar m, 
   vector Updf, 
   vector usc,
   scalar psi
   );

  //- Update end position (for tracking)
  void updateParticleProperties();

  //- Evolve the particles
  void evolve();

  //- Report Info
  void info() const;
  void particleInfo() const;
  void oneParticleInfo(const mcParticle& p) const;
  void assertPopulationHealth() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mcParticleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
