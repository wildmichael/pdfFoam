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

#include "mcSteadyFlamelet.H"

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
    chiIdx_(findIdx("chiName", "chi"))
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
                "mcSteadyFlamelet::mcSteadyFlamelet(mcParticleCloud&,"
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
                "rho",
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
    // sanity check: z_ strictly monothonically increasing
    for(label i = 1; i != Nz_; ++i)
    {
        if (!(z_[i]-z_[i-1] > 0))
        {
            FatalErrorIn("mcSteadyFlamelet::mcSteadyFlamelet("
                "mcParticleCloud&, const objectRegistry&, const word&)")
                << z_.objectPath()
                << " not strictly monothonically increasing.\n"
                << exit(FatalError);
        }
    }
    // sanity check: chi strictly monothonically increasing
    for(label i = 1; i != Nchi_; ++i)
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
    const volScalarField& zzCov =
        db().lookupObject<volScalarField>(zName_ + zName_ + "Cov");
    zVarInterp_ = interpolation<scalar>::New
    (
        cloud().solutionDict().interpolationScheme(zzCov.name()),
        zzCov
    );
}


void Foam::mcSteadyFlamelet::correct(mcParticle& p)
{
    const scalar& z = p.Phi()[zIdx_];
    const label& c = p.cell();
    const scalar zVar = zVarInterp_().interpolate(p.position(), c, p.face());
    const scalar chi = max(Cchi_*p.Omega()*zVar, SMALL);

    // compute index into z_ and interpolation weight
    label iz;
    scalar wz;
    findCoeffs(z_, z, iz, wz);
    // find index into chi_ and compute interpolation weight
    label ichi;
    scalar wchi;
    findCoeffs(chi_, chi, ichi, wchi);
    p.Phi()[chiIdx_] = chi;
    // interpolate rho
    p.rho() = binterp(phi_[0], ichi, wchi, iz, wz);
    // interpolate user-data
    forAll(addedIdx_, i)
    {
        const scalarListIOList& v = phi_[nNativeFields_+i];
        p.Phi()[addedIdx_[i]] =
            binterp(v, ichi, wchi, iz, wz);
    }
}

// ************************************************************************* //
