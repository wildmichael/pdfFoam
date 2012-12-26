/*---------------------------------------------------------------------------*\
                pdfFoam: General Purpose PDF Solution Algorithm
                   for Reactive Flow Simulations in OpenFOAM

 Copyright (C) 2012 Michael Wild, Heng Xiao, Patrick Jenny,
                    Institute of Fluid Dynamics, ETH Zurich
-------------------------------------------------------------------------------
License
    This file is part of pdfFoam.

    pdfFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) version 3 of the same License.

    pdfFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with pdfFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::mcParticleCloud::mesh() const
{
    return mesh_;
}


inline const Foam::dimensionedScalar& Foam::mcParticleCloud::deltaT() const
{
    return deltaT_;
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


inline Foam::scalar Foam::mcParticleCloud::openingAngle() const
{
    return openingAngle_;
}


inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::mcParticleCloud::volumeOrArea() const
{
    if (isAxiSymmetric_)
    {
        return area_();
    }
    else
    {
        return mesh_.V();
    }
}


inline const Foam::surfaceVectorField& Foam::mcParticleCloud::CourantCoeffs() const
{
    return CourantCoeffs_;
}


inline Foam::scalar
Foam::mcParticleCloud::massPerDepth(const mcParticle& p) const
{
    if (isAxiSymmetric_)
    {
        scalar r = mag(p.position() - (p.position()&axis_)*axis_);
        return p.m()/(max(r, SMALL)*openingAngle_);
    }
    else
    {
        return p.m();
    }
}


template<class TrackData>
void Foam::mcParticleCloud::hitPatch
(
    mcParticle&       p,
    TrackData&        td,
#if FOAM_HEX_VERSION < 0x200
    const label       patchI
#else
    const label       patchI,
    const scalar      trackFraction,
    const tetIndices& tetIs
#endif
)
{
#if FOAM_HEX_VERSION < 0x200
    boundaryHandlers_[patchI].hitPatch(p, td);
#else
    boundaryHandlers_[patchI].hitPatch(p, td, patchI, trackFraction, tetIs);
#endif
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
