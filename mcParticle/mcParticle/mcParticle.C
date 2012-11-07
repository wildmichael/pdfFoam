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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcParticle::trackData::trackData(mcParticleCloud& mcpc, scalar trackTime)
:
    Particle<mcParticle>::trackData(mcpc),
    cloud_(mcpc),
    trackTime_(trackTime),
    deltaT_(mcpc.deltaT().value())
{}


Foam::mcParticle::mcParticle
(
    const mcParticleCloud& c,
    const vector& position,
    const label   celli,
    const scalar  m,
    const vector& UParticle,
    const scalarField&  Phi,
    const vector& shift,
    const label   ghost
)
:
    Particle<mcParticle>(c, position, celli),
    m_(m),
    UParticle_(UParticle),
    Ucorrection_(vector::zero),
    Utracking_(UParticle),
    Omega_(0.0),
    rho_(0.0),
    eta_(1.0),
    shift_(shift),
    Co_(0.0),
    reflectionBoundaryVelocity_(vector::zero),
    ghost_(ghost),
    nSteps_(0),
    isOnInletBoundary_(false),
    reflectedAtOpenBoundary_(false),
    Phi_(Phi)
{
    const polyMesh& mesh = c.mesh();
    meshTools::constrainDirection(mesh, mesh.geometricD(), Utracking_);
    c.computeCourantNo(*this);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mcParticle::move(mcParticle::trackData& td)
{

    td.switchProcessor = false;
    td.keepParticle = true;

    mcParticleCloud& mcpc = td.cloud();

    if (isOnInletBoundary_)
    {
        stepFraction() = mcpc.random().scalar01();
    }

    const polyMesh& mesh = mcpc.pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    scalar trackTime = eta_*td.trackTime();
    scalar tEnd = (1.0 - stepFraction())*trackTime;
    scalar dtMax = tEnd;

    isOnInletBoundary_ = false;

    while (td.keepParticle && !td.switchProcessor && tEnd > 0)
    {
        if (debug)
        {
            Pout<< "Time = " << mesh.time().timeName()
                << "  trackTime = " << trackTime
                << "  tEnd = " << tEnd
                << "  stepFraction() = " << stepFraction()
                << "  origId() = " << origId()
                << "  position() = " << position();
        }

        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);
        point destPos = position() + dt * Utracking_;
        meshTools::constrainDirection(mesh, mesh.geometricD(), destPos);
        scalar tf = trackToFace(destPos, td);
        ++nSteps_;

        // if we made too many very small steps, drop the particle
        if (nSteps_ > 1000)
        {
#ifdef FULLDEBUG
            Perr<< "DEBUG: particle " << origId_ << " made more than 1000 "
                   "steps, droping it. Info:\n"
                << info() << nl;
#endif
            mcpc.notifyLostParticle(*this);
            td.keepParticle = false;
            break;
        };

        dt *= tf;
        if (debug)
        {
            Pout<< "  trackFraction = " << tf
                << endl;
        }

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/trackTime;

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
    Ucorrection_ = transform(T, Ucorrection_);
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
        << "Ucorrection = " << Ucorrection() << ", "
        << "Utracking = " << Utracking() << ", "
        << "U     = " << UParticle()  << ", "
        << "Phi   = " << Phi() << ", "
        << "ghost = " << ghost() << ", "
        << "shift = " << shift()
        << endl;
    return oss.str();
}


// ************************************************************************* //
