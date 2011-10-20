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
// #include "interpolationCellPointFaceFlux.H"
#include "boundBox.H"
#include "fvc.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineParticleTypeNameAndDebug(mcParticle, 0);
    defineTemplateTypeNameAndDebug(Cloud<mcParticle>, 0);
};


inline bool Foam::mcParticleCloud::less::operator()
(
 mcParticle* one,
 mcParticle* two
) const
{
  return one->m() < two->m();
}

inline bool Foam::mcParticleCloud::more::operator()
(
 mcParticle* one,
 mcParticle* two
) const
{
  return one->m() > two->m();
}


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
    debug_ (false),
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
    Npc_(particleProperties_.lookupOrAddDefault<label>("particlesPerCell", 20)),
    histNp_(size()),

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

    instantM0_(M0_*1.0),

    ghostCellHash_(256),
    ghostFaceHash_(256),

    UcPdf_
    (
      IOobject
     (
      "UCloudPDF",
      mesh.time().timeName(),
      mesh,
       IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
       ),
      M1_/max(M0_, SMALL_MASS)
     ),
    
    TaucPdf_
    (
      IOobject
      (
       "TauCloudPDF",
       mesh.time().timeName(),
       mesh,
       IOobject::READ_IF_PRESENT,
       IOobject::AUTO_WRITE
       ),
     M2_/max(M0_, SMALL_MASS)
     )

{

  if (size() > 0) // if particle data not found
    {
      mcParticle::readFields(*this);
    }
  else
    {
      Info << "I am releasing particles initially!" << endl;
      initReleaseParticles();
    }
  
  // Take care of statistical moments (make sure they are consistent)
  checkMoments();
  // Ensure particles takes the updated PDF values
  updateParticlePDF();

  findGhostLayers();

  mesh_.time().write();
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
  else if(size() > 0)
    {
      Info << "Moments are missing. Forced re-initialization." << endl;
      updateCloudPDF(0.0);
    }

}


// Perform particle averaging to obtained cell-based values.
void Foam::mcParticleCloud::updateCloudPDF(scalar existWt)
{
  instantM0_ *= 0.0;
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

      instantM0_[p.cell()] += p.m();
      instantM1[p.cell()] += p.m() * p.UParticle();
      instantM2[p.cell()] += p.m() * symm(u * u);
    }

  M0_ = existWt * M0_ + (1.0 - existWt) * instantM0_;
  M1_ = existWt * M1_ + (1.0 - existWt) * instantM1;
  M2_ = existWt * M2_ + (1.0 - existWt) * instantM2;
  // Compute U and tau
  UcPdf_  = M1_/max(M0_, SMALL_MASS);
  TaucPdf_  = M2_/max(M0_, SMALL_MASS);
}


//- update Updf & (maybe later, Taupdf) in the particle side.
void Foam::mcParticleCloud::updateParticlePDF()
{
  for(iterator pIter=begin(); 
      pIter != end();
      ++pIter
      )
    {
      mcParticle & p = pIter();
      p.Updf() = UcPdf_[p.cell()];
      // p.Taupdf() = TaucPdf_[p.cell()];
    }
}


// Update particle-cell addressing: for each cell, what are the IDs of the
// particles I host?
void Foam::mcParticleCloud::particleNumberControl()
{

  List<cellPopStatus> cellPopFlag(mesh_.nCells(), NORMAL);

  forAll(PaNIC_, celli)
    {
      label np = PaNIC_[celli];

      // classify the particle # heath condition of each cell
      if (np < 1)
        { cellPopFlag[celli] = EMPTY; }
      else if ( np <= (Npc_ * 2 / 3) )
        { cellPopFlag[celli] = TOOFEW; }
      else if (np >= Npc_* 3 / 2)
        { cellPopFlag[celli] = TOOMANY; }
      
      // clear old list
      cellParticleAddr_[celli].clear(); 
    }

    labelList ncpi(mesh_.nCells(), 0);

    for(mcParticleCloud::iterator pIter=begin(); 
      pIter != end();
      ++pIter
      )
    {
      mcParticle & p = pIter();
      label celli = p.cell();
      // If a mcParticleList is necessary for this cell, construct it.
      if(cellPopFlag[celli] > NORMAL) // either too few or too many particles
        { 
          // according to ascending order of mass (if too many particles)
          // or descending order of mass (if too few particle)
          mcParticleList & cepl = cellParticleAddr_[celli];
          if(cepl.size() < 1)  cepl.setSize(PaNIC_[celli]);
          cepl[ncpi[celli]++] = &p;
        }
    }
  
    // Sort the lists with particles, and perform particle number control
    forAll(cellPopFlag, celli)
      { 

        if (cellPopFlag[celli] == TOOFEW)
          {      
            // Exclude ghost cells
            if (!ghostCellHash_.found(celli))
              {
                sort(cellParticleAddr_[celli], more());
                cloneParticles(celli);
              }
          }
        else if ( cellPopFlag[celli] == TOOMANY )
          {
            if (!ghostCellHash_.found(celli))
              {
                sort(cellParticleAddr_[celli],  less());
                clusterParticles(celli);
              }
          }
      }

        
  // Debug only: check the list
    if(debug_)
      forAll(cellPopFlag, celli)
        {
          if (cellPopFlag[celli] == TOOFEW || cellPopFlag[celli] == TOOMANY )
            { 
              Info << "size is " << cellParticleAddr_[celli].size()
                   << ", flag is " << cellPopFlag[celli]
                   << " := " << endl;
              mcParticleList & cepl = cellParticleAddr_[celli];
              
              forAll(cepl, pci) 
                {
                  Info << "ID= " << cepl[pci]->origId() << ", " 
                       << "m= " << cepl[pci]->m() << endl;
                }
            }
        }
}

// Split the n heaviest particles
void Foam::mcParticleCloud::cloneParticles(label celli)
{
  label n = Npc_ - PaNIC_[celli]; // # particle to eliminate
  n = min(PaNIC_[celli], n);

  for(label particleI=0; particleI < n; particleI++)
    {
      mcParticle& p = *(cellParticleAddr_[celli][particleI]);
      // Half my mass
      p.m() /= 2.0;
      // create a new particle like myself
      autoPtr<mcParticle> ptrNew = p.clone();

      addParticle(ptrNew.ptr());
    }

  histNp_ += n;
}


// As name suggests
void Foam::mcParticleCloud::clusterParticles(label celli)
{
  label ncur = PaNIC_[celli];
  label nx =  ncur - Npc_; // # particle to eliminate
  // Pool of partiles to operate on is 2*nx, 
  // but liminted by available particles in this cell.  
  label nPool = min(2*nx+4, ncur); 
  
  label nKilled = 0;
  scalar massKilled = 0.0;
  for(label particleI=0; particleI < nPool; particleI++)
    {
      mcParticle& p = * cellParticleAddr_[celli][particleI];
      if ( random().scalar01() > 0.5 )
        { 
          massKilled += p.m();
          deleteParticle(p);
          cellParticleAddr_[celli][particleI] = 0;
          nKilled++;
        }
      if (nKilled >= nx) break;
    }

  scalar massCompensate = massKilled / (ncur - nKilled);

  // Compensate for the deleted mass
  for (label particleI=0; particleI < ncur; particleI++)
    {
      mcParticle* pPtr =  cellParticleAddr_[celli][particleI];
      if(pPtr)
        { 
          pPtr->m() += massCompensate; 
        }
    }

}


// Find "ghost cells" (actually the first layer of cells of in/out-flow patch
void Foam::mcParticleCloud::findGhostLayers()
{
  const cellList& cells = mesh_.cells();
  const faceList& faces = mesh_.faces();
  const polyBoundaryMesh& patches = mesh_.boundaryMesh();
  // in/out-flow patches (should read from dictionary)
  const wordList patchNames(particleProperties_.lookup("inoutPatches"));

  label ngPatchs = patchNames.size();
  labelList patchSizes(ngPatchs);

  if(ngPatchs > 0)
    {
      ghostCellLayers_.setSize(ngPatchs);
      ghostFaceLayers_.setSize(ngPatchs);
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
          FatalErrorIn("mcParticleCloud::findGhostLayers()")
            << "Illegal patch " << patchNames[nameI]
            << ". Valid patches are " << patches.name()
            << exit(FatalError);
        }
      label nf = patches[patchI].size();

      patchSizes[nameI] = nf;
      
      if (nf > 0)
        {
          ghostCellLayers_[nameI].setSize(nf);
          ghostFaceLayers_[nameI].setSize(nf);
          
          // Find the ghost cells and faces
          const polyPatch& curPatch = patches[patchI];
          label j = 0;
          label k = 0;
          forAll(curPatch, facei)
            {
              // find ghost cell
              label faceCelli = curPatch.faceCells()[facei];
              ghostCellLayers_[nameI][j++] =  faceCelli;
                //- Add  to hash set
              ghostCellHash_.insert(faceCelli);

              // find opposite face
              label gFaceI = curPatch.start() + facei;
              const cell& ownCell = cells[faceCelli];
              label oppositeFaceI = ownCell.opposingFaceLabel(gFaceI, faces);
              if (oppositeFaceI == -1)
                {
                  FatalErrorIn("mcParticleCloud::findGhostLayers()")
                    << "Face:" << facei << " owner cell:" << faceCelli
                    << " is not a hex?" << abort(FatalError);
                }
              else
                {
                  ghostFaceLayers_[nameI][k++] =  oppositeFaceI;
                }

                //- Add  to hash set
              ghostFaceHash_.insert(oppositeFaceI);
            }
        }
    }
  // Check faces found
  
  /* forAll(ghostFaceLayers_, patchI)
     forAll(ghostFaceLayers_[patchI], facei)
   {
      Info <<  "face #: " << ghostFaceLayers_[patchI][facei] << endl;
      Info << faces[ghostFaceLayers_[patchI][facei]].normal(mesh_.points()) << endl;
      Info << faces[ghostFaceLayers_[patchI][facei]].centre(mesh_.points()) << endl;
    } */

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcParticleCloud::evolve()
{
     dimensionedScalar coeffCorr("coeffCorr", dimLength*dimLength/dimTime, 1e-4);
    const volScalarField& rho = mesh_.lookupObject<const volScalarField>("rho");
    const volVectorField& U = mesh_.lookupObject<const volVectorField>("U");
    const volVectorField& gradP = mesh_.lookupObject<const volVectorField>("grad(p)");
    const volScalarField& k = mesh_.lookupObject<const volScalarField>("k");
    const volScalarField& epsilon = mesh_.lookupObject<const volScalarField>("epsilon");
    volScalarField diffRho(M0_);
    diffRho.internalField() = -(M0_.internalField()/mesh_.V() - rho)/rho;
    volVectorField gradRho = fvc::grad(diffRho) * coeffCorr;

    //    Info << "current time index: " << mesh_.time().timeIndex() << endl;

    interpolationCellPoint<scalar> rhoInterp(rho);
    //interpolationCellPointFaceFlux UInterp(U);
    interpolationCellPoint<vector> UInterp(U);
    interpolationCellPoint<vector> gradPInterp(gradP);
    interpolationCellPoint<scalar> kInterp(k);
    interpolationCellPoint<scalar> epsilonInterp(epsilon);
    interpolationCellPoint<vector> gradRhoInterp(gradRho);

    mcParticle::trackData td(*this, rhoInterp, UInterp, gradPInterp, kInterp, 
                             epsilonInterp, gradRho);

    Cloud<mcParticle>::move(td);

    // AvgTimeScale is a class member (should be read from dictionary).
    scalar existWt = 0.8;
      // 1.0/(1.0 + (mesh_.time().deltaT()/AvgTimeScale_).value());

    updateCloudPDF(existWt);
    updateParticlePDF();

    particleNumberControl();

    // Correct boundary conditions:
    populateGhostCells();

    //  particleInfo();
    assertPopulationHealth();
}


// Initialization: populate the FV field with particles
void Foam::mcParticleCloud::initReleaseParticles()
{
  // Populate each cell with partilces with N particle each cell
  label N = Npc_* 1.5;

  forAll(Ufv_, celli)
    {
      // &&& Should be a cloud property (class number)
      scalar m = mesh_.V()[celli] * rhofv_[celli] / N;
      vector Updf = Ufv_[celli];
      vector uscales(sqrt(kfv_[celli]), sqrt(kfv_[celli]), sqrt(kfv_[celli]));

      particleGenInCell(celli, N, m, Updf, uscales);
    }

  writeFields();
}


// Randomly generate N particles in celli, with provided cell-based
// values and the scale of velocity fluctuation
void Foam::mcParticleCloud::particleGenInCell
(
 label celli, 
 label N, 
 scalarList masses, 
 vector Updf, 
 vectorList uscales
 )
{
  boundBox cellbb(pointField(mesh_.points(), mesh_.cellPoints()[celli]));

  vector minb = cellbb.min();
  vector dimb = cellbb.max() - minb;

  label Npgen = 0;
  for(int i = 0; i < 100 * N; i++)
    {
      // Relative coordinate [0, 1] in this cell
      vector xi = random().vector01();
      // Random offset from min point
      scalar rx = min(max(10.0*SMALL, xi.x()), 1.0-10.0*SMALL);
      scalar ry = min(max(10.0*SMALL, xi.y()), 1.0-10.0*SMALL);
      scalar rz = min(max(10.0*SMALL, xi.z()), 1.0-10.0*SMALL);
      vector offsetRnd(rx*dimb.x(), ry*dimb.y(), rz*dimb.z());
      
      // Generate a particle position
      vector position = minb + offsetRnd;
      
      // Initially put $N particle per cell
      if(mesh_.pointInCell(position, celli))
        { 
          // What is the distribution of u?
          // How to enforce the component-wise correlations < u_i, u_j >?
          scalar m = masses[Npgen]; // &&& debug
          vector u( 
                    random().GaussNormal() * uscales[Npgen].x(),  
                    random().GaussNormal() * uscales[Npgen].y(),  
                    random().GaussNormal() * uscales[Npgen].z()
                  );
          vector UParticle = u + Updf;
          vector UFap = Ufv_[celli];
          mcParticle* ptr =
            new mcParticle
            (
             *this,  position, celli, m, Updf, UParticle, UFap, dtCloud_
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
  else // update count
    {
      histNp_ += Npgen;
    }
  
}


// Generate a number of particle with the same properties 
// (particularly: m and uscales), celli, Updf are always the same.
void Foam::mcParticleCloud::particleGenInCell
(
 label celli, 
 label N, 
 scalar m, 
 vector Updf, 
 vector usc
 )
{
  scalarList masses(N, m);
  vectorList uscales(N, usc);
  particleGenInCell(celli, N, masses, Updf, uscales);
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
              scalar m = (mesh_.V()[celli] * rhofv_[celli] - instantM0_[celli])  / N;
              if (m <= 0) 
                {
                  Info << "populateGhostCells::warning: negative mass"
                       << "patch " << patchi << "; "
                       << "cell " << celli << endl;
                  continue; // prevent negative mass
                }
              vector Updf = Ufv_[celli];
              scalar ksqrt = sqrt(kfv_[celli]);
              vector uscales(ksqrt, ksqrt, ksqrt);

              particleGenInCell(celli, N, m, Updf, uscales);
              if (debug_)
                Info << N << " particles generated in cell " << celli
                     << " m= " << m
                     << " original total mass:" << M0_[celli] 
                     << " oritinal avg mass: " <<  M0_[celli] / PaNIC_[celli]
                     << " total mass now = " << m * N + M0_[celli] << endl;
              
            }
        }
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
      << "    Current number of particles         = "
      << returnReduce(this->size(), sumOp<label>()) << nl
      << "    Totally number of particles existed = "
      << returnReduce(histNp_,      sumOp<label>()) << nl;
}


void Foam::mcParticleCloud::particleInfo() const
{
  for(mcParticleCloud::const_iterator pIter=begin(); 
      pIter != end();
      ++pIter
      )
    {
      const mcParticle & p = pIter();
      oneParticleInfo(p);
    }
}


void Foam::mcParticleCloud::oneParticleInfo(const mcParticle& p) const
{
  Pout << "Particle Id: " << p.origId() << ": "
       << "X    = " << p.position() << ", "
       << "cell = " << p.cell() << ", "
       << "m    = " << p.m() << nl
       << "U    = " << p.UParticle()  << ", "
       << "Ufv  = " << Ufv_[p.cell()] << ", "
       << "Updf = " << p.Updf()
       << endl;
}


void Foam::mcParticleCloud::assertPopulationHealth() const
{
  bool cloudHealthy = true;
  for(mcParticleCloud::const_iterator pIter=begin(); 
      pIter != end();
      ++pIter
      )
    {
      const mcParticle & p = pIter();
      if(! mesh_.bounds().contains(p.position()))
        {
          cloudHealthy = false;
          Info << "***Warning: " << endl;
          oneParticleInfo(p);
        }
      
    }

  if(!cloudHealthy)
    {
      FatalErrorIn("mcParticleCloud::assertPopulationHealth()") 
        << "Particle health check not passed!"
        << exit(FatalError);
    }

}
// ************************************************************************* //
