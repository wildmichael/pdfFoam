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

\*---------------------------------------------------------------------------*/

#include "mcParticleCloud.H"
#include "fvMesh.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "boundBox.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineParticleTypeNameAndDebug(mcParticle, 0);
    defineTemplateTypeNameAndDebug(Cloud<mcParticle>, 0);
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcParticleCloud::mcParticleCloud
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& k,
    const word& cloudName,
    bool readFields
)
:
    Cloud<mcParticle>(mesh, cloudName, false),
    mesh_(mesh),
    Ufv_(U),
    rhofv_(rho),
    kfv_(k),

    particleProperties_
    (
        IOobject
        (
            "particleProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
     ),

    dtCloud_(mesh.time().deltaT().value()),
    AvgTimeScale_(mesh.time().endTime().value()),
    random_(55555+12345*Pstream::myProcNo()),
    Npc_(10),

    M0_(
        IOobject
        (
            "M0",
            mesh_.time().constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
         ),
        mesh_,
        dimensionedScalar("M0", dimMass, 0.0)
       ),
    M1_(
        IOobject
        (
            "M1",
            mesh_.time().constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
         ),
        mesh_,
        dimensionedVector("M1", dimMass*dimVelocity, vector::zero)
        ),
    M2_(
        IOobject
          (
            "M2",
            mesh_.time().constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
           ),
        mesh_,
        dimensionedSymmTensor("M2", dimEnergy, symmTensor:: zero)
        ),

    UcPdf_(M1_/max(M0_, SMALL)),
    TaucPdf_(M2_/max(M0_, SMALL))
{
  //&&& should be a property of cloud (read from dict)
  bool initialRelease = true;
  bool readParticles = false; 

  if(initialRelease)
    {
      initReleaseParticles();
      readParticles = false;
    }
  
  if (readParticles)
    {
      mcParticle::readFields(*this);
    }
 
  // Take care of statistical moments (make sure they are consistent)
  checkMoments();

  // Ensure particles takes the updated PDF values
  updateParticlePDF();
 
}


//- update Updf & Taupdf in the particle side.
void Foam::mcParticleCloud::updateParticlePDF()
{
  for(iterator pIter=begin(); 
      pIter != end();
      ++pIter
      )
    {
      mcParticle & p = pIter();
      p.Updf() = UcPdf_[p.cell()];
      p.Taupdf() = TaucPdf_[p.cell()];
    }
}


// If moments are not read correctly, initialize them.
void Foam::mcParticleCloud::checkMoments()
{
  if( M0_.headerOk() &&
      M1_.headerOk() &&
      M2_.headerOk() )
    {
      Info << "Moments read correctly." << endl;
    }
  else
    {
      Info << "Moments are missing. Forced re-initialization." << endl;
      updateCloudPDF(0.0);
    }

}


// Perform particle averaging to obtained cell-based values.
void Foam::mcParticleCloud::updateCloudPDF(scalar existWt)
{
  volScalarField instantM0 = M0_ * 0.0;
  volVectorField instantM1 = M1_ * 0.0;
  volSymmTensorField instantM2 = M2_ * 0.0;

      // Loop through particles to accumulate moments (0, 1, 2 order)
      for(iterator pIter=begin(); 
          pIter != end();
          ++pIter
          )
        {
          mcParticle & p = pIter();
          vector u = p.UParticle() - p.Updf();
          instantM0[p.cell()] += p.m();
          instantM1[p.cell()] += p.m() * p.UParticle();
          instantM2[p.cell()] += p.m() * symm(u * u);
        }

      M0_ = existWt * M0_ + (1.0 - existWt) * instantM0;
      M1_ = existWt * M1_ + (1.0 - existWt) * instantM1;
      M2_ = existWt * M2_ + (1.0 - existWt) * instantM2;
      // Compute U and tau
      UcPdf_  = M1_/max(M0_, SMALL);
      TaucPdf_  = M2_/max(M0_, SMALL);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcParticleCloud::evolve()
{
    const volScalarField& rho = mesh_.lookupObject<const volScalarField>("rho");
    const volVectorField& U = mesh_.lookupObject<const volVectorField>("U");
    const volVectorField& gradP = mesh_.lookupObject<const volVectorField>("grad(p)");
    const volScalarField& k = mesh_.lookupObject<const volScalarField>("k");
    const volScalarField& epsilon = mesh_.lookupObject<const volScalarField>("epsilon");

    interpolationCellPoint<scalar> rhoInterp(rho);
    interpolationCellPoint<vector> UInterp(U);
    interpolationCellPoint<vector> gradPInterp(gradP);
    interpolationCellPoint<scalar> kInterp(k);
    interpolationCellPoint<scalar> epsilonInterp(epsilon);

    mcParticle::trackData td(*this, rhoInterp, UInterp, gradPInterp, kInterp, epsilonInterp);

    Cloud<mcParticle>::move(td);

    // AvgTimeScale is a class member (should be read from dictionary).
    scalar existWt = 
      1.0/(1.0 + (mesh_.time().deltaT()/AvgTimeScale_).value());
    
    updateCloudPDF(existWt);
    updateParticlePDF();

    allParticlesInfo();
}


// Initialization: populate the FV field with particles
void Foam::mcParticleCloud::initReleaseParticles()
{
  // Populate each cell with partilces
  forAll(Ufv_, celli)
    {
      // &&& Should be a cloud property (class number)
      scalar m = mesh_.V()[celli] * rhofv_[celli] / Npc_;
      vector Updf = UcPdf_[celli];
      vector uscales(sqrt(kfv_[celli]), sqrt(kfv_[celli]), sqrt(kfv_[celli]));

      particleGenInCell(celli, Npc_, m, Updf, uscales);
    }

  writeFields();
}


// Randomly generate N particles in celli, with provided cell-based
// values and the scale of velocity fluctuation
void Foam::mcParticleCloud::particleGenInCell
(
 label celli, 
 label N, 
 scalar m, 
 vector Updf, 
 vector uscales
 )
{
      // Generate $Npc_  particles randomly in celli
      // should be put in a function particleGenInCell(celli, Npc_)
      boundBox cellbb(pointField(mesh_.points(), mesh_.cellPoints()[celli]));

      vector minb = cellbb.min();
      vector dimb = cellbb.max()-minb;

      label Npgen = 0;
      for(int i = 0; i < 100 * Npc_; i++)
        {
          // Relative coordinate [0, 1] in this cell
          vector xi = random().vector01();
          // Random offset from min
          vector offsetRnd(xi.x()*dimb.x(), xi.y()*dimb.y(), xi.z()*dimb.z());
          
          // Generate a particle position
          vector position = minb + offsetRnd;
          
          /* Info << "cell #" << celli <<"; bb=  " << cellbb 
               << "particle generated: " << position
               << " In cell? " << mesh_.pointInCell(position, celli)  
               << endl; */
          
          // Initially put $Npc_ particle per cell
          if(mesh_.pointInCell(position, celli))
            { 
              // What is the distribution of u?
              // How to enforce the component-wise correlations < u_i, u_j >?
              vector u( random().GaussNormal() * uscales.x(),  
                        random().GaussNormal() * uscales.y(),  
                        random().GaussNormal() * uscales.z()
                       );
              vector UParticle = u + Updf;
              mcParticle* ptr =
                new mcParticle 
                (
                 *this,  position, celli, m, Updf, UParticle, dtCloud_
                 );

              addParticle(ptr);
              Npgen ++;
            }
          if(Npgen >= Npc_) break;
        }
      
      if(Npgen < Npc_) 
        {
          FatalErrorIn("mcParticleCloud::initReleaseParticles()")
            << "Only " << Npgen << " particles generated for cell "
            << celli << nl << "Something is wrong" << exit(FatalError);
        }
  
}


//This serves as a template for looping through particles in the cloud
void Foam::mcParticleCloud::updateParticleProperties()
{
 
  for(mcParticleCloud::iterator pIter=begin(); 
      pIter != end();
      ++pIter
      )
    {
      // mcParticle & p = pIter();
    }
  
}

void Foam::mcParticleCloud::info() const
{
  Info<< "Cloud: " << this->name() << nl
      << "    Current number of particles       = "
      << returnReduce(this->size(), sumOp<label>()) << nl;
}


void Foam::mcParticleCloud::allParticlesInfo() const
{
  Info << "Calling allParticlesInfo" << endl;
  for(mcParticleCloud::const_iterator pIter=begin(); 
      pIter != end();
      ++pIter
      )
    {
      const mcParticle & p = pIter();
      Pout << "Particle # " << p.origId() << "; "
           << "X = " << p.position() << "; "
           << "U = " << p.UParticle() << nl
           << "   Ufv = " << Ufv_[p.cell()] << "; "
           << "Updf = " << p.Updf()
           << endl;
      
    }
}

// ************************************************************************* //
