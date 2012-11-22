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

#include "mcKulkarniAutoIgnition.H"

#include "addToRunTimeSelectionTable.H"
#include "mcParticleCloud.H"
#include "interpolation.H"
#include "Pair.H"
#include "uniqueOrder_FIX.H"

#include <algorithm>

// * * * * * * * * * * * * * Local Helper Functions  * * * * * * * * * * * * //

namespace // anonymous
{

using namespace Foam;

//- Interpolation index and weight
template<class L>
inline void findCoeffs(const L& lst, scalar v, int& i, scalar& w)
{
    typename L::const_iterator iter =
        std::lower_bound(lst.cbegin(), lst.cend(), v);
    if (iter == lst.cbegin())
    {
        i = 0;
        w = 0.;
    }
    else if (iter == lst.cend())
    {
        i = lst.size() - 2;
        w = 1.;
    }
    else
    {
        typename L::const_iterator prev = iter;
        std::advance(prev, -1);
        i = std::distance(lst.cbegin(), prev);
        w = (v - *prev)/(*iter - *prev);
    }
    if (w < -VSMALL || w > (1.+VSMALL))
    {
        FatalErrorIn("::<anonymous>::findCoeffs<L>(...)")
            << "Weight w = " << w << " not within [0, 1]\n"
            << abort(FatalError);
    }
}

//- bilinear interpolation
template<class V>
inline scalar binterp
(
    const V& v, //!< data field (list-of-list, row-major)
    label ix,   //!< row index
    scalar wx,  //!< row weight
    label iy,   //!< first column index
    scalar wy   //!< column weight
)
{
    const scalar
        v00 = v[ix  ][iy  ],
        v10 = v[ix  ][iy+1],
        v01 = v[ix+1][iy  ],
        v11 = v[ix+1][iy+1];
    const scalar wxm  = 1. - wx;
    const scalar wym  = 1. - wy;
    return v00*wxm*wym + v10*wx*wym + v01*wxm*wy + v11*wx*wy;
}

} // anonymous namespace

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcKulkarniAutoIgnitionReactionModel, 0);
    addNamedToRunTimeSelectionTable
    (
        mcReactionModel,
        mcKulkarniAutoIgnitionReactionModel,
        mcReactionModel,
        KulkarniAutoIgnition
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcKulkarniAutoIgnitionReactionModel::mcKulkarniAutoIgnitionReactionModel
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcReactionModel(cloud, db, subDictName),
    zName_  (thermoDict().lookupOrDefault<word>("zName", "z")),
    cName_  (thermoDict().lookupOrDefault<word>("cName", "z")),
    z_
    (
        IOobject
        (
            zName_,
            db.time().constant(),
            "autoIgnition",
            db,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    pv_
    (
        IOobject
        (
            "pv",
            db.time().constant(),
            "autoIgnition",
            db,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    cEq_
    (
        IOobject
        (
            "cEq",
            db.time().constant(),
            "autoIgnition",
            db,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    Nz_(z_.size()),
    Npv_(pv_.size()),
    cEqMax_(max(cEq_)),
    addedNames_(),
    phi_(),
    zIdx_(findIdx("zName", "z")),
    cIdx_(findIdx("cName", "c")),
    pvIdx_(findIdx("pvName", "pv"))
{
    if (thermoDict().found("scalars"))
    {
        thermoDict().lookup("scalars") >> addedNames_;
        labelList order;
        uniqueOrder_FIX(addedNames_, order);
        if (addedNames_.size() != order.size())
        {
            FatalErrorIn
            (
                "mcKulkarniAutoIgnitionReactionModel::"
                "mcKulkarniAutoIgnitionReactionModel(mcParticleCloud&,"
                "const objectRegistry&, const word&)"
            )
                << "The list "
                << thermoDict().lookupEntry("scalars", false, false).name()
                << " contains duplicate entries.\n"
                << exit(FatalError);
        }
        addedIdx_.setSize(addedNames_.size());
        forAll(addedNames_, i)
        {
            word n = addedNames_[i];
            addedIdx_[i] = findIdx(n+"Name", n);
        }
    }
    phi_.setSize(nNativeFields_ + addedNames_.size());
    // load density constant data
    phi_.set
    (
        0,
        new scalarListIOList
        (
            IOobject
            (
                "cdot",
                db.time().constant(),
                "autoIgnition",
                db,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        )
    );
    phi_.set
    (
        1,
        new scalarListIOList
        (
            IOobject
            (
                "rho",
                db.time().constant(),
                "autoIgnition",
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
                    "autoIgnition",
                    db,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );
    }
    // z and pv must have at least two elements
    if (Nz_ < 2)
    {
        FatalErrorIn("mcKulkarniAutoIgnitionReactionModel::"
            "mcKulkarniAutoIgnitionReactionModel("
            "mcParticleCloud&, const objectRegistry&, const word&)")
            << z_.objectPath() << " must have at least 2 elements.\n"
            << exit(FatalError);
    }
    if (Npv_ < 2)
    {
        FatalErrorIn("mcKulkarniAutoIgnitionReactionModel::"
            "mcKulkarniAutoIgnitionReactionModel("
            "mcParticleCloud&, const objectRegistry&, const word&)")
            << pv_.objectPath() << " must have at least 2 elements.\n"
            << exit(FatalError);
    }
    // cEq must have the same length as z
    if (cEq_.size() != Nz_)
    {
        FatalErrorIn("mcKulkarniAutoIgnitionReactionModel::"
            "mcKulkarniAutoIgnitionReactionModel("
            "mcParticleCloud&, const objectRegistry&, const word&)")
            << cEq_.objectPath() << " must have the same size as "
            << z_.objectPath() << ".\n"
            << exit(FatalError);
    }
    // sanity check: z_ strictly monothonically increasing
    for(label i = 1; i != Nz_; ++i)
    {
        if (!(z_[i]-z_[i-1] > 0))
        {
            FatalErrorIn("mcKulkarniAutoIgnitionReactionModel::mcKulkarniAutoIgnitionReactionModel("
                "mcParticleCloud&, const objectRegistry&, const word&)")
                << z_.objectPath()
                << " not strictly monothonically increasing.\n"
                << exit(FatalError);
        }
    }
    // sanity check: pv strictly monothonically increasing
    for(label i = 1; i != Npv_; ++i)
    {
        if (!(pv_[i]-pv_[i-1] > 0))
        {
            FatalErrorIn("mcKulkarniAutoIgnitionReactionModel::mcKulkarniAutoIgnitionReactionModel("
                "mcParticleCloud&, const objectRegistry&, const word&)")
                << pv_.objectPath()
                << " not strictly monothonically increasing.\n"
                << exit(FatalError);
        }
    }
    // sanity check: phi[*].size() == Nz_ and phi[*][*].size() == Npv_
    forAll(phi_, i)
    {
        if (phi_[i].size() != Nz_)
        {
            FatalErrorIn("mcKulkarniAutoIgnitionReactionModel::"
                "mcKulkarniAutoIgnitionReactionModel("
                "mcParticleCloud&, const objectRegistry&, const word&)")
                << phi_[i].objectPath()
                << " must have the same number of entries as "
                << z_.objectPath() << " (" << Nz_ << ")\n"
                << exit(FatalError);
        }
        forAll(phi_[i], j)
        {
            if (phi_[i][j].size() != Npv_)
            {
                FatalErrorIn("mcKulkarniAutoIgnitionReactionModel::"
                    "mcKulkarniAutoIgnitionReactionModel("
                    "mcParticleCloud&, const objectRegistry&, const word&)")
                    << phi_[i].objectPath()
                    << " entry " << i << " must have the same number of"
                    << " entries as " << pv_.objectPath()
                    << " (" << Npv_ << ")\n"
                    << exit(FatalError);
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::mcKulkarniAutoIgnitionReactionModel::correct(mcParticle& p)
{
    const scalar& deltaT = cloud().deltaT().value();
    const scalar& z = p.Phi()[zIdx_];
    scalar& c = p.Phi()[cIdx_];

    // compute index into z_ and interpolation weight
    label iz;
    scalar wz;
    findCoeffs(z_, z, iz, wz);

    // equilibirum progress variable
    scalar cEq = cEq_[iz]*(1. - wz) + cEq_[iz+1]*wz;

    // limit c to cEq
    c = min(c, cEq);

    // compute normalised progress variable
    scalar pv = cEq > 1e-5*cEqMax_ ? c/cEq : 0.;
    p.Phi()[pvIdx_] = pv;

    // find index into pv_ and compute interpolation weight
    label ic;
    scalar wc;
    findCoeffs(pv_, pv, ic, wc);

    // interpolate cdot and integrate in time
    scalar cdot = binterp(phi_[0], iz, wz, ic, wc);
    c += cdot*p.eta()*deltaT;

    // update weights and interpolate rho
    findCoeffs(pv_, pv, ic, wc);
    p.rho() = binterp(phi_[1], iz, wz, ic, wc);

    // interpolate user-data
    forAll(addedIdx_, i)
    {
        p.Phi()[addedIdx_[i]] =
            binterp(phi_[nNativeFields_+i], iz, wz, ic, wc);
    }
    p.Co() = max(p.Co(), cdot);
}

// ************************************************************************* //
