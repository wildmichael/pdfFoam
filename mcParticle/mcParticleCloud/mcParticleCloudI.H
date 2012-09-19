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

inline const Foam::fvMesh& Foam::mcParticleCloud::mesh() const
{
    return mesh_;
}


inline const Foam::dictionary& Foam::mcParticleCloud::thermoDict() const
{
    return thermoDict_;
}


inline const Foam::mcSolution& Foam::mcParticleCloud::solutionDict() const
{
    return solutionDict_;
}


inline bool Foam::mcParticleCloud::read()
{
    solutionDict_.read();
    return regIOobject::read();
}


inline bool Foam::mcParticleCloud::readIfModified()
{
    solutionDict_.readIfModified();
    return regIOobject::readIfModified();
}


inline const Foam::compressible::turbulenceModel&
Foam::mcParticleCloud::turbulenceModel() const
{
    return turbModel_;
}


inline const Foam::labelList& Foam::mcParticleCloud::mixedScalars() const
{
    return mixedScalars_;
}


inline const Foam::labelList& Foam::mcParticleCloud::conservedScalars() const
{
    return conservedScalars_;
}


inline const Foam::volVectorField& Foam::mcParticleCloud::Ufv() const
{
    return Ufv_;
}

inline const Foam::tmp<Foam::volScalarField> Foam::mcParticleCloud::kfv() const
{
    return turbModel_.k();
}


inline const Foam::tmp<Foam::volScalarField>
Foam::mcParticleCloud::epsilonfv() const
{
    return turbModel_.epsilon();
}

inline const Foam::volScalarField& Foam::mcParticleCloud::rhocPdfInst() const
{
    return rhocPdfInst_;
}

inline const Foam::volScalarField& Foam::mcParticleCloud::rhocPdf() const
{
    return rhocPdf_;
}

inline const Foam::volScalarField& Foam::mcParticleCloud::kcPdf() const
{
    return kcPdf_;
}


inline const Foam::List<Foam::volScalarField*>&
Foam::mcParticleCloud::PhicPdf() const
{
    return PhicPdf_;
}


inline const Foam::volSymmTensorField&
Foam::mcParticleCloud::TaucPdf() const
{
    return TaucPdf_;
}


inline const Foam::volScalarField&
Foam::mcParticleCloud::pndcPdf() const
{
    return pndcPdf_;
}


inline const Foam::volScalarField&
Foam::mcParticleCloud::pndcPdfInst() const
{
    return pndcPdfInst_;
}


inline const Foam::wordList& Foam::mcParticleCloud::scalarNames() const
{
    return scalarNames_;
}


inline bool Foam::mcParticleCloud::isAxiSymmetric() const
{
    return isAxiSymmetric_;
}


inline const Foam::vector& Foam::mcParticleCloud::axis() const
{
    return axis_;
}


inline const Foam::vector& Foam::mcParticleCloud::centrePlaneNormal() const
{
    return centrePlaneNormal_;
}


inline const Foam::surfaceVectorField& Foam::mcParticleCloud::CourantCoeffs() const
{
    return CourantCoeffs_;
}


template<class TrackData>
void Foam::mcParticleCloud::hitPatch
(
    Foam::mcParticle& p,
    TrackData& td,
    Foam::label patchI
)
{
    boundaryHandlers_[patchI].hitPatch(p, td);
}


inline const Foam::mcLocalTimeStepping&
Foam::mcParticleCloud::localTimeStepping() const
{
    return localTimeStepping_();
}


inline Foam::mcLocalTimeStepping&
Foam::mcParticleCloud::localTimeStepping()
{
    return localTimeStepping_();
}

// ************************************************************************* //
