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


inline const Foam::tmp<Foam::volScalarField> Foam::mcParticleCloud::kfv() const
{
    return turbModel_.k();
}


inline const Foam::tmp<Foam::volScalarField>
Foam::mcParticleCloud::epsilonfv() const
{
    return turbModel_.epsilon();
}


inline const Foam::volScalarField& Foam::mcParticleCloud::kcPdf() const
{
    return kcPdf_;
}

inline const Foam::dimensionedScalar& Foam::mcParticleCloud::kRelaxTime() const
{
    return kRelaxTime_;
}

inline Foam::label Foam::mcParticleCloud::Npc()
{
    return Npc_;
}

// ************************************************************************* //
