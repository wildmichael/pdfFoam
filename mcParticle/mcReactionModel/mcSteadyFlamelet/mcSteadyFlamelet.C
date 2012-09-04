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

#include "mcSteadyFlamelet.H"

#include "addToRunTimeSelectionTable.H"
#include "mcParticleCloud.H"
#include "Pair.H"

#include <algorithm>

// * * * * * * * * * * * * * Local Helper Functions  * * * * * * * * * * * * //

namespace // anonymous
{

using namespace Foam;

//- bilinear interpolation
template<class V>
inline scalar binterp
(
    const V& v, //!< data field (list-of-list, row-major)
    label ix1,  //!< first row index
    label ix2,  //!< second row index
    scalar wx,  //!< row weight
    label iy1,  //!< first column index
    label iy2,  //!< second column index
    scalar wy   //!< column weight
)
{
    const scalar
        v00 = v[ix1][iy1],
        v10 = v[ix1][iy2],
        v01 = v[ix2][iy1],
        v11 = v[ix2][iy2];
    const scalar wxm  = 1. - wx;
    const scalar wym  = 1. - wy;
    return v00*wxm*wym + v10*wx*wym + v01*wxm*wy + v11*wx*wy;
}

} // anonymous namespace

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcSteadyFlamelet, 0);
    addNamedToRunTimeSelectionTable
    (
        mcReactionModel,
        mcSteadyFlamelet,
        mcReactionModel,
        SteadyFlamelet
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcSteadyFlamelet::mcSteadyFlamelet
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcReactionModel(cloud, db, subDictName),
    zName_  (thermoDict().lookupOrDefault<word>("zName", "z")),
    TName_(thermoDict().lookupOrDefault<word>("TName", "T")),
    Cchi_(thermoDict().lookupOrDefault<scalar>("Cchi", 6.0)),
    z_
    (
        IOobject
        (
            zName_,
            db.time().constant(),
            "flamelet",
            db,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    chi_
    (
        IOobject
        (
            "chi",
            db.time().constant(),
            "flamelet",
            db,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    Nz_(z_.size()),
    Nchi_(chi_.size()),
    addedNames_(),
    phi_(),
    zIdx_(findIdx("zName", "z")),
    chiIdx_(findIdx("chiName", "chi")),
    TIdx_(findIdx("TName", "T")),
    p_
    (
        db.lookupObject<volScalarField>
            (thermoDict().lookupOrDefault<word>("pName", "p"))
    )
{
    if (thermoDict().found("scalars"))
    {
        thermoDict().lookup("scalars") >> addedNames_;
        addedIdx_.setSize(addedNames_.size());
        forAll(addedNames_, i)
        {
            word n = addedNames_[i];
            addedIdx_[i] = findIdx(n+"Name", n);
        }
    }
    phi_.setSize(nNativeFields_ + addedNames_.size());
    // load temperature data
    phi_.set
    (
        0,
        new scalarListIOList
        (
            IOobject
            (
                TName_,
                db.time().constant(),
                "flamelet",
                db,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        )
    );
    // load specific gas constant data
    phi_.set
    (
        1,
        new scalarListIOList
        (
            IOobject
            (
                "Rgas",
                db.time().constant(),
                "flamelet",
                db,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        )
    );
    // load user data
    forAll(addedNames_, i)
    {
        phi_.set
        (
            nNativeFields_ + i,
            new scalarListIOList
            (
                IOobject
                (
                    addedNames_[i],
                    db.time().constant(),
                    "flamelet",
                    db,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );
    }
    // z must have at least two, chi one element
    if (Nz_ < 2)
    {
        FatalErrorIn("mcSteadyFlamelet::mcSteadyFlamelet("
            "mcParticleCloud&, const objectRegistry&, const word&)")
            << z_.objectPath() << " must have at least 2 elements.\n"
            << exit(FatalError);
    }
    if (Nchi_ < 1)
    {
        FatalErrorIn("mcSteadyFlamelet::mcSteadyFlamelet("
            "mcParticleCloud&, const objectRegistry&, const word&)")
            << chi_.objectPath() << " must have at least 1 element.\n"
            << exit(FatalError);
    }
    // sanity checks: deltaZ ~ constant, z_ strictly monothonically increasing,
    // at the same time compute mean deltaZ
    scalar dz = z_[1] - z_[0];
    bool incr = dz > 0;
    deltaZ_ = dz;
    for(label i = 2; i != Nz_ && incr; ++i)
    {
        scalar dzn = z_[i] - z_[i-1];
        if (abs(dzn-dz)/dz > 1e-3)
        {
            FatalErrorIn("mcSteadyFlamelet::mcSteadyFlamelet("
                "mcParticleCloud&, const objectRegistry&, const word&)")
                << z_.objectPath() << " not uniformly spaced (up to 0.1%).\n"
                << exit(FatalError);
        }
        dz = dzn;
        deltaZ_ += dz;
        incr &= dz > 0;
    }
    deltaZ_ /= Nz_;
    if (!incr)
    {
        FatalErrorIn("mcSteadyFlamelet::mcSteadyFlamelet("
            "mcParticleCloud&, const objectRegistry&, const word&)")
            << z_.objectPath() << " not strictly monothonically increasing.\n"
            << exit(FatalError);
    }
    // sanity check: chi strictly monothonically increasing
    for(label i = 1; i != Nchi_ && incr; ++i)
    {
        if (!(chi_[i]-chi_[i-1] > 0))
        {
            FatalErrorIn("mcSteadyFlamelet::mcSteadyFlamelet("
                "mcParticleCloud&, const objectRegistry&, const word&)")
                << chi_.objectPath()
                << " not strictly monothonically increasing.\n"
                << exit(FatalError);
        }
    }
    // sanity check: phi[*].size() == Nchi_ and phi[*][*].size() == Nz_
    forAll(phi_, i)
    {
        if (phi_[i].size() != Nchi_)
        {
            FatalErrorIn("mcSteadyFlamelet::mcSteadyFlamelet("
                "mcParticleCloud&, const objectRegistry&, const word&)")
                << phi_[i].objectPath()
                << " must have the same number of entries as "
                << chi_.objectPath() << " (" << Nchi_ << ")\n"
                << exit(FatalError);
        }
        forAll(phi_[i], j)
        {
            if (phi_[i][j].size() != Nz_)
            {
                FatalErrorIn("mcSteadyFlamelet::mcSteadyFlamelet("
                    "mcParticleCloud&, const objectRegistry&, const word&)")
                    << phi_[i].objectPath()
                    << " entry " << i << " must have the same number of"
                    << " entries as " << z_.objectPath()
                    << " (" << Nz_ << ")\n"
                    << exit(FatalError);
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcSteadyFlamelet::updateInternals()
{
    mcReactionModel::updateInternals();
    zVarInterp_.reset
    (
        new interpolationCellPointFace<scalar>
        (
            db().lookupObject<volScalarField>(zName_ + zName_ + "Cov")
        )
    );
}


void Foam::mcSteadyFlamelet::correct(mcParticle& p)
{
    const scalar& z = p.Phi()[zIdx_];
    const label& c = p.cell();
    const scalar zVar = zVarInterp_().interpolate(p.position(), c, p.face());
    const scalar chi = max(Cchi_*p.Omega()*zVar, SMALL);

    // compute index into z_ and interpolation weight
    label iz1 = floor(z/deltaZ_);
    label iz2;
    scalar wz = 0;
    if (iz1 < 0)
    {
        iz1 = iz2 = 0;
    }
    else if (iz1 > Nz_ - 2)
    {
        iz1 = iz2 = Nz_ - 1;
    }
    else
    {
        iz2 = iz1 + 1;
        wz = (z - z_[iz1])/(z_[iz2] - z_[iz1]);
    }
    // find index into chi_ and compute interpolation weight
    const scalarList::const_iterator chiIter =
        std::lower_bound(chi_.cbegin(), chi_.cend(), chi);
    label ichi1, ichi2;
    scalar wchi = 0;
    if (chiIter == chi_.cbegin())
    {
        ichi1 = ichi2 = 0;
    }
    else if (chiIter > chi_.cend() - 2)
    {
        ichi1 = ichi2 = Nchi_ - 1;
    }
    else
    {
        ichi1 = std::distance(chi_.cbegin(), chiIter);
        ichi2 = ichi1 + 1;
        wchi = (chi - chi_[ichi1])/(chi_[ichi2] - chi_[ichi1]);
    }
    // interpolate T and R
    const scalar T = binterp(phi_[0], ichi1, ichi2, wchi, iz1, iz2, wz);
    const scalar R = binterp(phi_[1], ichi1, ichi2, wchi, iz1, iz2, wz);
    p.Phi()[chiIdx_] = chi;
    p.Phi()[TIdx_] = T;
    // interpolate user-data
    forAll(addedIdx_, i)
    {
        const scalarListIOList& v = phi_[nNativeFields_+i];
        p.Phi()[addedIdx_[i]] =
            binterp(v, ichi1, ichi2, wchi, iz1, iz2, wz);
    }
    // compute particle density
    p.rho() = p_[c]/(R*T);
}

// ************************************************************************* //
