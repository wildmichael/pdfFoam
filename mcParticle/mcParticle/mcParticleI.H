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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolationCellPointFace<Foam::scalar>&
Foam::mcParticle::trackData::rhoInterp() const
{
    return rhoInterp_;
}


inline const Foam::interpolationCellPointFace<Foam::vector>&
Foam::mcParticle::trackData::UInterp() const
{
    return UInterp_;
}


inline const Foam::gradInterpolationConstantTet<Foam::scalar>&
Foam::mcParticle::trackData::gradPInterp() const
{
    return gradPInterp_;
}


inline const Foam::interpolationCellPointFace<Foam::scalar>&
Foam::mcParticle::trackData::kInterp() const
{
    return kInterp_;
}


inline const Foam::interpolationCellPointFace<Foam::vector>&
Foam::mcParticle::trackData::diffUInterp() const
{
    return diffUInterp_;
}


inline Foam::scalar Foam::mcParticle::m() const
{
    return m_;
}


inline Foam::scalar& Foam::mcParticle::m()
{
    return m_;
}


inline const Foam::vector& Foam::mcParticle::Updf() const
{
    return Updf_;
}


inline Foam::vector& Foam::mcParticle::Updf()
{
    return Updf_;
}


inline const Foam::vector& Foam::mcParticle::UParticle() const
{
    return UParticle_;
}


inline Foam::vector& Foam::mcParticle::UParticle()
{
    return UParticle_;
}


inline const Foam::vector& Foam::mcParticle::Ucorrection() const
{
    return Ucorrection_;
}


inline Foam::vector& Foam::mcParticle::Ucorrection()
{
    return Ucorrection_;
}


inline const Foam::vector& Foam::mcParticle::Utracking() const
{
    return Utracking_;
}


inline const Foam::vector& Foam::mcParticle::UFap() const
{
    return UFap_;
}


inline Foam::scalar Foam::mcParticle::Omega() const
{
  return Omega_;
}


inline Foam::scalar& Foam::mcParticle::Omega()
{
  return Omega_;
}


inline const Foam::scalarField& Foam::mcParticle::Phi() const
{
  return Phi_;
}


inline Foam::scalarField& Foam::mcParticle::Phi()
{
  return Phi_;
}


inline Foam::scalar Foam::mcParticle::rho() const
{
  return rho_;
}


inline Foam::scalar& Foam::mcParticle::rho()
{
  return rho_;
}

inline Foam::scalar Foam::mcParticle::eta() const
{
  return eta_;
}

inline Foam::scalar& Foam::mcParticle::eta()
{
  return eta_;
}

inline Foam::label Foam::mcParticle::ghost() const
{
    return ghost_;
}


inline Foam::label& Foam::mcParticle::ghost()
{
    return ghost_;
}


inline Foam::label Foam::mcParticle::nSteps() const
{
    return nSteps_;
}


inline Foam::label& Foam::mcParticle::nSteps()
{
    return nSteps_;
}


inline const Foam::vector& Foam::mcParticle::shift() const
{
    return shift_;
}


inline Foam::vector& Foam::mcParticle::shift()
{
    return shift_;
}


inline Foam::scalar Foam::mcParticle::Co() const
{
    return Co_;
}


inline const Foam::vector& Foam::mcParticle::reflectionBoundaryVelocity() const
{
    return reflectionBoundaryVelocity_;
}


inline Foam::vector& Foam::mcParticle::reflectionBoundaryVelocity()
{
    return reflectionBoundaryVelocity_;
}


inline bool Foam::mcParticle::isOnInletBoundary() const
{
    return isOnInletBoundary_;
}


inline bool& Foam::mcParticle::isOnInletBoundary()
{
    return isOnInletBoundary_;
}


inline bool Foam::mcParticle::reflectedAtOpenBoundary() const
{
    return reflectedAtOpenBoundary_;
}


inline bool& Foam::mcParticle::reflectedAtOpenBoundary()
{
    return reflectedAtOpenBoundary_;
}

// ************************************************************************* //
