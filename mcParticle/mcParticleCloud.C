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
    )
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
    const volScalarField& k = mesh_.lookupObject<const volScalarField>("k");

    interpolationCellPoint<scalar> rhoInterp(rho);
    interpolationCellPoint<vector> UInterp(U);
    interpolationCellPoint<scalar> kInterp(k);

    mcParticle::trackData td(*this, rhoInterp, UInterp, kInterp);

    updateEndPositions();

    Cloud<mcParticle>::move(td);
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

      if (celli < 800) continue;
      // Initially put $nci particle per cell
      for(int i = 0; i < nci_; i++)
        {
          vector u = vector(-Updf.x()+1, 0, 0);
          
          mcParticle* ptr = new mcParticle(*this,  position, celli, m, Updf, u, dtCloud_);
          
          addParticle(ptr);
        }

    }
  
}


void Foam::mcParticleCloud::updateEndPositions()
{
 
  for(mcParticleCloud::iterator pIter=begin(); 
      pIter != end();
      ++pIter
      )
    {
      mcParticle & p = pIter();
      p.endPosition() = p.position() + (p.UParticle()) * p.dt();
    }
  
}

void Foam::mcParticleCloud::info() const
{
  Info<< "Cloud: " << this->name() << nl
      << "    Current number of particles       = "
      << returnReduce(this->size(), sumOp<label>()) << nl;
}



// ************************************************************************* //
