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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::mcParticle::trackData::trackData
(
    mcParticleCloud& mcpc,
    const interpolationCellPoint<scalar>& rhoInterp,
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCellPoint<scalar>& kInterp
)
:
    mcpc_(mcpc),
    rhoInterp_(rhoInterp),
    UInterp_(UInterp),
    kInterp_(kInterp)
{}


inline Foam::mcParticle::mcParticle
(
 const Cloud<mcParticle>& c,
 const vector& position,
 const label celli,
 const scalar m,
 const vector& Updf,
 const vector& u,
 const scalar dt
)
:
  Particle<mcParticle>(c, position, celli),
  m_(m),
  Updf_(Updf),
  u_(u),
  endPosition_(position),
  dt_(dt)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::mcParticleCloud& Foam::mcParticle::trackData::mcpc()
{
    return mcpc_;
}

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::mcParticle::trackData::rhoInterp() const
{
    return rhoInterp_;
}

inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::mcParticle::trackData::UInterp() const
{
    return UInterp_;
}

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::mcParticle::trackData::kInterp() const
{
    return kInterp_;
}


inline Foam::scalar Foam::mcParticle::m() const
{
    return m_;
}

inline const Foam::vector& Foam::mcParticle::Updf() const
{
    return Updf_;
}

inline const Foam::vector& Foam::mcParticle::u() const
{
    return u_;
}

inline const Foam::vector Foam::mcParticle::UParticle()
{
    return Updf_ + u_;
}

inline Foam::scalar Foam::mcParticle::dt()
{
    return dt_;
}

inline Foam::vector & Foam::mcParticle::endPosition()
{
  return endPosition_;
}

// ************************************************************************* //
