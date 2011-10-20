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
    const word& cloudName,
    bool readFields
)
:
    Cloud<mcParticle>(mesh, cloudName, false),
    mesh_(mesh),
    Ufv_(U),
    dtCloud_(mesh.time().deltaT().value()),
    random_(55555+12345*Pstream::myProcNo()),
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
    M0_(
        IOobject
        (
            "M0",
            mesh_.time().constant(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
        mesh_,
        
       ),
    M1_(),
    M2_()
{
  //&&& should be a property cloud (read from dict)
  bool initialRelease = true;

  if(initialRelease)
    {
      initReleaseParticles();
      readFields = false;
    }
  
  if (readFields)
    {
      mcParticle::readFields(*this);
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
    allParticlesInfo();
}


// Initialization: populate the FV field with particles
void Foam::mcParticleCloud::initReleaseParticles()
{
  // Populate each cell with partilces
  forAll(Ufv_, celli)
    {
      // &&& Should be a cloud property.
      label nci_ = 1;
      scalar m = mesh_.V()[celli]/nci_;
      vector position = mesh_.C()[celli];
      vector Updf = Ufv_[celli];

      if (celli !=  400) continue;
      // Initially put $nci particle per cell
      for(int i = 0; i < nci_; i++)
        {
          vector u = vector(2, 1 , 0);
          
          mcParticle* ptr = new mcParticle(*this,  position, celli, m, Updf, u, dtCloud_);
          
          addParticle(ptr);
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
      Info << "Particle # " << p.origId() << "; "
           << "X = " << p.position() << "; "
           << "Up = " << p.UParticle()
           << endl;
      
    }
}

// ************************************************************************* //
