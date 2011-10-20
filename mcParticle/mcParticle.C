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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mcParticle::move(mcParticle::trackData& td)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    while (td.keepParticle && !td.switchProcessor 
           && stepFraction_ < 1.0 - SMALL )
    {
        if (debug)
        {
            Info<< "Time = " << mesh.time().timeName()
                << " steptFraction() = " << stepFraction() << endl;
        }

        stepFraction_ += trackToFace(endPosition(), td)*(1.0 - stepFraction_);  

        if (onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[patch(face())]))
              {
                td.switchProcessor = true;
              }
            //else if (isA<wallPolyPatch>(pbMesh[patch(face())]))
            // {
                // Info << "Hitting wall patch..." << endl;
                //const wallPolyPatch & wpp = static_cast<const wallPolyPatch&> (pbMesh[patch(face())]);
                //hitWallPatch(wpp, td);
            // }
        }
    }

    return td.keepParticle;
}


bool Foam::mcParticle::hitPatch
(
 const polyPatch& patch,
 mcParticle::trackData& td,
 const label Lb
)
{
  bool skipPatchHitAction = false;
  
  if(
     isA<wallPolyPatch>(patch) ||
     isA<processorPolyPatch>(patch)
     )
    {
      Info << "Known patch type" << endl;
    }
  else
    {
      Info << "Unknown patch type" << endl;
    }
    
  return skipPatchHitAction;
}


bool Foam::mcParticle::hitPatch
(
    const polyPatch&,
    int&,
    const label
)
{
  return false;
}


void Foam::mcParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    mcParticle::trackData& td
)
{
    td.switchProcessor = true;
}


void Foam::mcParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    int&
)
{}


void Foam::mcParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    mcParticle::trackData& td
)
{
  vector nw = wpp.faceAreas()[wpp.whichFace(face())];
  nw /= mag(nw);  // Wall normal (outward)
  
  vector UParticle = Updf_ + u_;
  
  scalar Un = UParticle & nw; // Normal component
  vector Ut = UParticle - Un*nw; // Tangential component
  
  if (Un > 0)
    {
      // Elastic coefficient: 1 = perfectly elastic; 0 = perfect plastic.
      UParticle -= (1.0 + td.mcpc().e())*Un*nw; 
    }

  // Tangential friction coefficient: 
  // 1 = totally friction (stop); 0 = smooth (no loss)
  UParticle -= td.mcpc().mu()*Ut;
  
  u_ = UParticle - Updf_;  // Updf is not changed.
}


void Foam::mcParticle::hitWallPatch
(
    const wallPolyPatch&,
    int&
)
{}


void Foam::mcParticle::hitPatch
(
    const polyPatch&,
    mcParticle::trackData& td
)
{
  Info << "action for unknown patch type." << endl;
  td.keepParticle = false;
}


void Foam::mcParticle::hitPatch
(
    const polyPatch&,
    int&
)
{}


void Foam::mcParticle::transformProperties (const tensor& T)
{
    Particle<mcParticle>::transformProperties(T);
    vector UParticle = Updf_ + u_;
    UParticle = transform(T, UParticle);
    u_ = UParticle - Updf_;
}


void Foam::mcParticle::transformProperties(const vector& separation)
{
    Particle<mcParticle>::transformProperties(separation);
}


// ************************************************************************* //
