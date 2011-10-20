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
            IOobject::AUTO_WRITE
        )
     ),

    dtCloud_(mesh.time().deltaT().value()),
    AvgTimeScale_(mesh.time().endTime().value()),
    random_(55555+12345*Pstream::myProcNo()),
    Npc_(50),

    SMALL_MASS("SMALL_MASS", dimMass, SMALL),

    PaNIC_(
        IOobject
        (
         "PaNIC",
         mesh.time().timeName(),
         mesh_,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
         ),
        mesh.nCells()
        ),

    cellParticleAddr_(mesh.nCells()),

    M0_(
        IOobject
        (
            "M0",
            mesh.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
         ),
        mesh,
        dimensionedScalar("M0", dimMass, 0.0)
       ),
    M1_(
        IOobject
        (
            "M1",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
         ),
        mesh,
        dimensionedVector("M1", dimMass*dimVelocity, vector::zero)
        ),
    M2_(
        IOobject
          (
            "M2",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
           ),
        mesh,
        dimensionedSymmTensor("M2", dimEnergy, symmTensor:: zero)
        ),

    UcPdf_
    (
     // IOobject
     // (
     //  "UcloudPDF",
     //  mesh.time().timeName(),
     //  mesh,
     //  IOobject::READ_IF_PRESENT,
     //  IOobject::AUTO_WRITE
     //  ),
     // mesh,
     M1_/max(M0_, SMALL_MASS)
     ),
    
    TaucPdf_
    (
     // IOobject
     // (
     //  "TauCloudPDF",
     //  mesh.time().timeName(),
     //  mesh,
     //  IOobject::READ_IF_PRESENT,
     //  IOobject::AUTO_WRITE
     //  ),
     // mesh,
     M2_/max(M0_, SMALL_MASS)
     )

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

  updateCellParticleAddr();

  UcPdf_.writeOpt() = IOobject::AUTO_WRITE; 
  UcPdf_.rename("UcPdf");
  TaucPdf_.writeOpt() = IOobject::AUTO_WRITE; 
  TaucPdf_.rename("TauPdf");
  mesh_.time().write();

  findGhostCellLayers();

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
  if( PaNIC_.headerOk() &&
      M0_.headerOk() &&
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

  PaNIC_ *= 0;

  // Loop through particles to accumulate moments (0, 1, 2 order)
  // as well as particle number
  for(iterator pIter=begin(); 
      pIter != end();
      ++pIter
      )
    {
      mcParticle & p = pIter();
      vector u = p.UParticle() - p.Updf();
      PaNIC_[p.cell()] ++;
      instantM0[p.cell()] += p.m();
      instantM1[p.cell()] += p.m() * p.UParticle();
      instantM2[p.cell()] += p.m() * symm(u * u);
    }
  
  M0_ = existWt * M0_ + (1.0 - existWt) * instantM0;
  M1_ = existWt * M1_ + (1.0 - existWt) * instantM1;
  M2_ = existWt * M2_ + (1.0 - existWt) * instantM2;
  // Compute U and tau
  UcPdf_  = M1_/max(M0_, SMALL_MASS);
  TaucPdf_  = M2_/max(M0_, SMALL_MASS);
}


// Update particle-cell addressing: for each cell, what are the IDs of the
// particles I host?
void Foam::mcParticleCloud::updateCellParticleAddr()
{

  labelList npi(mesh_.nCells(), 0);

  forAll(PaNIC_, celli)
    {
      label np = PaNIC_[celli];
      if (np > 0)
        cellParticleAddr_[celli].setSize(np); 
      else
        cellParticleAddr_[celli].clear(); 
    }

  for(iterator pIter=begin(); 
      pIter != end();
      ++pIter
      )
    {
      mcParticle & p = pIter();
      label celli = p.cell();
      // NOTES: it is not possible to use this origId to quickly locate
      // particles in the cloud (for deletion). This part needs more work. 
      cellParticleAddr_[npi[celli]] = p.origId(); 
      npi[celli] ++;
    }
  
}


// Find "ghost cells" (actually the first layer of cells of in/out-flow patch
void Foam::mcParticleCloud::findGhostCellLayers()
{
  const polyBoundaryMesh& patches = mesh_.boundaryMesh();
  
  // in/out-flow patches (should read from dictionary)
  const wordList patchNames(particleProperties_.lookup("inoutPatches"));

  label ngPatchs = patchNames.size();
  labelList patchSizes(ngPatchs);

  if(ngPatchs > 0)
    {
      ghostCellLayers_.setSize(ngPatchs);
    }
  else
    { 
      return; 
    }


  forAll(patchNames, nameI)
    {
      label patchI = patches.findPatchID(patchNames[nameI]);
      
      if (patchI == -1)
        {
          FatalErrorIn("mcParticleCloud::findGhostCellLayers()")
            << "Illegal patch " << patchNames[nameI]
            << ". Valid patches are " << patches.name()
            << exit(FatalError);
        }
      label nf = patches[patchI].size();

      patchSizes[nameI] = nf;
      
      if (nf > 0)
        {
          ghostCellLayers_[nameI].setSize(nf);
          
          // Find the ghost cells
          const polyPatch& curPatch = patches[patchI];
          label j = 0;
          forAll(curPatch, facei)
            {
              label faceCelli = curPatch.faceCells()[facei];
              ghostCellLayers_[nameI][j++] =  faceCelli;
            }
        }
    }
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
    scalar existWt = 0.8;
      // 1.0/(1.0 + (mesh_.time().deltaT()/AvgTimeScale_).value());
    
    updateCloudPDF(existWt);
    updateParticlePDF();
    updateCellParticleAddr();

    // Correct boundary conditions:
    populateGhostCells();

    //  allParticlesInfo();
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
  boundBox cellbb(pointField(mesh_.points(), mesh_.cellPoints()[celli]));

  vector minb = cellbb.min();
  vector dimb = cellbb.max()-minb;

  label Npgen = 0;
  for(int i = 0; i < 100 * N; i++)
    {
      // Relative coordinate [0, 1] in this cell
      vector xi = random().vector01();
      // Random offset from min point
      vector offsetRnd(xi.x()*dimb.x(), xi.y()*dimb.y(), xi.z()*dimb.z());
      
      // Generate a particle position
      vector position = minb + offsetRnd;
      
      // Initially put $N particle per cell
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
      
      // until enough particles are generated.
      if(Npgen >= N) break;
    }
  

  if(Npgen < N) 
    {
      FatalErrorIn("mcParticleCloud::initReleaseParticles()")
        << "Only " << Npgen << " particles generated for cell "
        << celli << nl << "Something is wrong" << exit(FatalError);
    }
  
}


// Enforce in/out flow BCs by populating ghost cells
void Foam::mcParticleCloud::populateGhostCells()
{
  forAll(ghostCellLayers_, patchi)
    {
      forAll(ghostCellLayers_[patchi], faceCelli)
        {
          label celli = ghostCellLayers_[patchi][faceCelli];
          
          if(PaNIC_[celli] < Npc_ / 2)
            { 
              label N = Npc_- PaNIC_[celli];
              scalar m = (mesh_.V()[celli] * rhofv_[celli] - M0_[celli])  / N;
              if (m <= 0) continue;
              vector Updf = Ufv_[celli];
              scalar ksqrt = sqrt(kfv_[celli]);
              vector uscales(ksqrt, ksqrt, ksqrt);
              //NOTE: the mass is not correct.
              particleGenInCell(celli,  N, m, Updf, uscales);
              Info << N << " particles generated in cell " << celli
                   << " m= " << m
                   << " original total mass:" << M0_[celli] 
                   << " oritinal avg mass: " <<  M0_[celli] / PaNIC_[celli]
                   << " total mass now = " << m * N + M0_[celli] << endl;
              
            }
        }
    }

  label globalNp = size();
  reduce(globalNp, sumOp<label>());
  Info << " Current particle number in the system: " << size()
       << "sum of all processors: " << globalNp << endl;
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
