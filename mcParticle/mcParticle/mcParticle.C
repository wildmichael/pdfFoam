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

#include "mcParticle.H"

#include "mcParticleCloud.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * Local Helper Functions  * * * * * * * * * * * * //

namespace // anonymous
{

Foam::scalar computeCourantNo(const Foam::mcParticle& p)
{
    using namespace Foam;
    const mcParticleCloud& cloud = refCast<const mcParticleCloud>(p.cloud());
    const polyMesh& mesh = cloud.pMesh();
    const surfaceVectorField& CourantCoeffs = cloud.CourantCoeffs();
    const cell& c = mesh.cells()[p.cell()];
    const vector& U = p.Utracking();
    scalar dt = mesh.time().deltaT().value();
    scalar Co = 0.0;
    forAll(c, cellFaceI)
    {
        label faceI = c[cellFaceI];
        vector coeff;
        if (mesh.isInternalFace(faceI))
        {
            coeff = CourantCoeffs[faceI];
        }
        else
        {
            label patchI = mesh.boundaryMesh().whichPatch(faceI);
            label start = mesh.boundaryMesh()[patchI].start();
            coeff = CourantCoeffs.boundaryField()[patchI][faceI-start];
        }
        Co = max(Co, fabs(dt * coeff & U));
    }
    return Co;
}

} // anonymous namespace

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcParticle::trackData::trackData
(
    mcParticleCloud& mcpc,
    const interpolationCellPoint<scalar>& rhoInterp,
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCellPoint<vector>& gradPInterp,
    const interpolationCellPoint<scalar>& kInterp,
    const interpolationCellPoint<vector>& gradRhoInterp,
    const interpolationCellPoint<vector>& diffUInterp
)
:
    Particle<mcParticle>::trackData(mcpc),
    rhoInterp_(rhoInterp),
    UInterp_(UInterp),
    gradPInterp_(gradPInterp),
    kInterp_(kInterp),
    gradRhoInterp_(gradRhoInterp),
    diffUInterp_(diffUInterp)
{}


Foam::mcParticle::mcParticle
(
    const mcParticleCloud& c,
    const vector& position,
    const label   celli,
    const scalar  m,
    const vector& Updf,
    const vector& UParticle,
    const vector& UFap,
    const scalarField&  Phi,
    const scalar  dt,
    const vector& shift,
    const label   ghost
)
:
    Particle<mcParticle>(c, position, celli),
    m_(m),
    Updf_(Updf),
    UParticle_(UParticle),
    UFap_(UFap),
    Omega_(0.0),
    rho_(0.0),
    dt_(dt),
    shift_(shift),
    Co_(0.0),
    ghost_(ghost),
    Phi_(Phi)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mcParticle::move(mcParticle::trackData& td)
{
    // SLM constant, temporarily put here C0 = 2.1
    const scalar C0 = 2.1;

    td.switchProcessor = false;
    td.keepParticle = true;

    mcParticleCloud& mcpc = refCast<mcParticleCloud>(td.cloud());

    const polyMesh& mesh = mcpc.pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    scalar deltaT = mesh.time().deltaT().value();
    scalar tEnd = (1.0 - stepFraction())*deltaT;
    scalar dtMax = tEnd;

    if (stepFraction() < SMALL)
    {
        if (debug)
        {
            Info<< "Updating fluctuating velocity" << endl;
        }
        cellPointWeight cpw(mesh, position(), cell(), face());

        // fluid quantities @ particle position
        vector gradPFap = td.gradPInterp().interpolate(cpw);
        scalar kFap = td.kInterp().interpolate(cpw);
        vector diffUap = td.diffUInterp().interpolate(cpw);

        // interpolate fluid velocity to particle location This
        // quantity is a data member the class.
        // Note: if would be the best if UInterp is interpolating
        // velocities based on faces to get UFap instead of cell
        // center values. Will implemente later.
        UFap_ = td.UInterp().interpolate(cpw);

        // Wiener process (question mark)
        vector dW = sqrt(deltaT) * vector
        (
            mcpc.random().GaussNormal(),
            mcpc.random().GaussNormal(),
            mcpc.random().GaussNormal()
        );

        // Update velocity
        UParticle_ += - gradPFap/rho_ * deltaT
            - (0.5 + 0.75 * C0) * Omega_ * (UParticle_- Updf_) * deltaT
            + sqrt(C0 * kFap * Omega_) * dW
            + diffUap * deltaT;

        // Scale to ensure consistency on TKE
        UParticle_ +=
              (UParticle_ - Updf_)  * deltaT / mcpc.kRelaxTime().value()
            * (sqrt(mcpc.kfv()()[cell()]/mcpc.kcPdf()[cell()]) - 1.0);
    }

    vector correctedUp = UParticle_ - Updf_ + UFap_;

    cellPointWeight cpwx(mesh, position(), cell(), face());
    vector gradRhoFap = td.gradRhoInterp().interpolate(cpwx);

    Utracking_ = correctedUp - gradRhoFap;
    meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking_);

    point destPos = position() + tEnd * Utracking_;

    if (mcpc.isAxiSymmetric())
    {
        vector rotatedCentreNormal = mcpc.axis() ^ destPos;
        rotatedCentreNormal /= mag(rotatedCentreNormal);
        tensor T = rotationTensor(rotatedCentreNormal, mcpc.centrePlaneNormal());
        transformProperties(T);
        destPos = transform(T, destPos);
        // constrain to kill numerical artifacts
        meshTools::constrainDirection(mesh, mesh.geometricD(), Utracking_);
        meshTools::constrainDirection(mesh, mesh.geometricD(), destPos);
    }

    Co_ = computeCourantNo(*this);

    while (td.keepParticle && !td.switchProcessor && tEnd > SMALL)
    {
        if (debug)
        {
            Info<< "Time = " << mesh.time().timeName()
                << "  deltaT = " << deltaT
                << "  tEnd = " << tEnd
                << "  stepFraction() = " << stepFraction()
                << "  origId() = " << origId()
                << "  position() = " << position();
        }

        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);
        destPos = position() + dt * Utracking_;

        // Particle does not move with its actual velocity, but with
        // FV interpolated velocity plus particle fluctuation
        // velocity. A particle number correction flux is needed as
        // well (to ensure consistency with FV density field).

        scalar tf = trackToFace(destPos, td);
        dt *= tf;
        if (debug)
        {
            Info<< "  trackFraction = " << tf
                << endl;
        }

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/deltaT;

        if (onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[patch(face())]))
            {
                td.switchProcessor = true;
            }
        }
    }

    return td.keepParticle;
}



// Pre-action before hitting patches
bool Foam::mcParticle::hitPatch
(
    const polyPatch& patch,
    mcParticle::trackData& td,
    const label Lb
)
{
    return false;
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

    scalar Un = UParticle_ & nw; // Normal component
    vector Ut = UParticle_ - Un*nw; // Tangential component

    if (Un > 0)
    {
        UParticle_ -= 2.0*Un*nw;
    }
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
    // Only transform fluctuating velocity
    UParticle_ = transform(T, UParticle_-Updf_)+Updf_;
    Utracking_ = transform(T, Utracking_);
}


void Foam::mcParticle::transformProperties(const vector& separation)
{
    Particle<mcParticle>::transformProperties(separation);
}


Foam::string Foam::mcParticle::info() const
{
    OStringStream oss;
    oss << "Particle Id: " << origId() << ": "
        << "X     = " << position() << ", "
        << "cell  = " << cell() << ", "
        << "m     = " << m() << nl
        << "U     = " << UParticle()  << ", "
        << "UFap  = " << UFap() << ", "
        << "Updf  = " << Updf() << nl
        << "Phi   = " << Phi() << ", "
        << "ghost = " << ghost() << ", "
        << "shift = " << shift()
        << endl;
    return oss.str();
}


// ************************************************************************* //
