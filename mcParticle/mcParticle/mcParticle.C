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
            if (!CourantCoeffs.boundaryField()[patchI].size())
            {
                // skip empty boundaries
                continue;
            }
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
    const interpolationCellPointFace<scalar>& rhoInterp,
    const interpolationCellPointFace<vector>& UInterp,
    const interpolationCellPointFace<vector>& gradPInterp,
    const interpolationCellPointFace<scalar>& kInterp,
    const interpolationCellPointFace<vector>& gradRhoInterp,
    const interpolationCellPointFace<vector>& diffUInterp
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
    shift_(shift),
    Co_(0.0),
    reflectionBoundaryVelocity_(vector::zero),
    ghost_(ghost),
    nSteps_(0),
    isOnInletBoundary_(false),
    reflectedAtOpenBoundary_(false),
    Phi_(Phi)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mcParticle::move(mcParticle::trackData& td)
{

    td.switchProcessor = false;
    td.keepParticle = true;

    mcParticleCloud& mcpc = refCast<mcParticleCloud>(td.cloud());
    const scalar& C0 = mcpc.C0();
    const scalar& C1 = mcpc.C1();

    if (isOnInletBoundary_)
    {
        stepFraction() = mcpc.random().scalar01();
    }

    const polyMesh& mesh = mcpc.pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    scalar deltaT = mesh.time().deltaT().value();
    scalar tEnd = (1.0 - stepFraction())*deltaT;
    scalar dtMax = tEnd;

    if (stepFraction() < SMALL || isOnInletBoundary_)
    {
        if (debug)
        {
            Info<< "Updating fluctuating velocity" << endl;
        }
        const point& p = position();
        label c = cell();
        label f = face();

        // fluid quantities @ particle position
        vector gradPFap = td.gradPInterp().interpolate(p, c, f);
        scalar kFap = td.kInterp().interpolate(p, c, f);
        vector diffUap = td.diffUInterp().interpolate(p, c, f);

        // interpolate fluid velocity to particle location This
        // quantity is a data member the class.
        // Note: if would be the best if UInterp is interpolating
        // velocities based on faces to get UFap instead of cell
        // center values. Will implemente later.
        UFap_ = td.UInterp().interpolate(p, c, f);

        // Wiener process (question mark)
        vector dW = sqrt(deltaT) * vector
        (
            mcpc.random().GaussNormal(),
            mcpc.random().GaussNormal(),
            mcpc.random().GaussNormal()
        );

        // Update velocity
        UParticle_ += - gradPFap/rho_ * deltaT
            - (0.5 * C1 + 0.75 * C0) * Omega_ * (UParticle_- UFap_) * deltaT
            + sqrt(C0 * kFap * Omega_) * dW
            // Correct mean velocity
            + diffUap * deltaT
            // Scale to ensure consistency on TKE (using interpolated
            // kFap/kPdfap gives very bad results)
            + (UParticle_ - UFap_)  * deltaT / mcpc.kRelaxTime().value()
              * (sqrt(mcpc.kfv()()[cell()]/mcpc.kcPdf()[cell()]) - 1.0);
    }

    isOnInletBoundary_ = false;

    vector correctedUp = UParticle_; // - Updf_ + UFap_;

    vector gradRhoFap = td.gradRhoInterp().interpolate
    (
        position(),
        cell(),
        face()
    );

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
        ++nSteps_;

        // if we made too many very small steps, drop the particle
        if (nSteps_ > 1000)
        {
#ifdef FULLDEBUG
            Perr<< "DEBUG: particle " << origId_ << " made more than 100 "
                   "steps, droping it. Info:\n"
                << info() << nl;
#endif
            mcpc.notifyLostParticle(*this);
            td.keepParticle = 0;
            break;
        };

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
#define DEFINE_HITPATCH(TRACKDATA)                                            \
bool Foam::mcParticle::hitPatch                                               \
(                                                                             \
    const polyPatch& patch,                                                   \
    TRACKDATA&       td,                                                      \
    const label      patchI                                                   \
)                                                                             \
{                                                                             \
    if (isA<wedgePolyPatch>(patch))                                           \
    {                                                                         \
        const polyMesh& mesh = cloud().pMesh();                               \
        meshTools::constrainDirection(mesh, mesh.geometricD(), position());   \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        mcParticleCloud& c = const_cast<mcParticleCloud&>                     \
        (                                                                     \
            refCast<const mcParticleCloud>(cloud())                           \
        );                                                                    \
        c.hitPatch(*this, td, patchI);                                        \
    }                                                                         \
    return true;                                                              \
}

DEFINE_HITPATCH(mcParticle::trackData)
DEFINE_HITPATCH(int)


void Foam::mcParticle::transformProperties (const tensor& T)
{
    Particle<mcParticle>::transformProperties(T);
    // Only transform fluctuating velocity
    UParticle_ = transform(T, UParticle_);
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
        << "Utracking = " << Utracking() << ", "
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
