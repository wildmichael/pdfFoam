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

#include "mcRASOmegaModel.H"

#include "addToRunTimeSelectionTable.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "compressible/RAS/kOmegaSST/kOmegaSST.H"
#include "interpolation.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * Local Helper Functions  * * * * * * * * * * * * //

namespace // anonymous
{


//- Similar to fvc::bound(), but clip to maximum instead of minimum
void boundMax(Foam::volScalarField& vsf, const Foam::dimensionedScalar vsf0)
{
    using namespace Foam;
    scalar maxVsf = gMax(vsf);
    if (maxVsf > vsf0.value())
    {
        Info<< "bounding " << vsf.name()
            << ", min = " << gMin(vsf) << ", max = " << maxVsf
            << ", average = " << gAverage(vsf.internalField()) << nl;
        vsf = min(vsf, vsf0);
    }
}


} // anonymous namespace

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcRASOmegaModel, 0);
    addNamedToRunTimeSelectionTable
    (
        mcOmegaModel,
        mcRASOmegaModel,
        mcOmegaModel,
        RAS
    );

} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcRASOmegaModel::mcRASOmegaModel
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcOmegaModel(cloud, db, subDictName),
    Omega_(),
    OmegaInterp_(),
    Omega0_("Omega0", dimless/dimTime, HUGE)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcRASOmegaModel::updateInternals()
{
    const compressible::turbulenceModel& turbModel = cloud().turbulenceModel();
    if (isA<compressible::RASModels::kOmegaSST>(turbModel))
    {
        // If we have a kOmegaSST (or derived) object, use omega directly
        const compressible::RASModels::kOmegaSST& kOmegaModel =
            refCast<const compressible::RASModels::kOmegaSST>(turbModel);
        Omega_.reset
        (
            new volScalarField("mcRASOmegaModel::Omega", kOmegaModel.omega())
        );
    }
    else
    {
        // Otherwise compute Omega from epsilon/k
        Omega_.reset
        (
            new volScalarField
            (
                "mcRASOmegaModel::Omega",
                turbModel.epsilon()/turbModel.k()
            )
        );
    }
    Omega0_.readIfPresent(solutionDict());
    boundMax(Omega_(), Omega0_);
    OmegaInterp_ = interpolation<scalar>::New
    (
        cloud().solutionDict().interpolationScheme(Omega_().name()),
        Omega_()
    );
}


void Foam::mcRASOmegaModel::correct(Foam::mcParticle& p)
{
    if (!Omega_.valid())
    {
        FatalErrorIn("mcRASOmegaModel::correct"
            "(mcParticleCloud&,mcParticle&,bool)")
            << "autoPtr holding Omega not valid"
            << exit(FatalError);
    }
    if (!OmegaInterp_.valid())
    {
        FatalErrorIn("mcRASOmegaModel::correct"
            "(mcParticleCloud&,mcParticle&,bool)")
            << "Interpolator not initialized"
            << exit(FatalError);
    }
    p.Omega() =
        max
        (
            OmegaInterp_().interpolate(p.position(), p.cell(), p.face()),
            SMALL
        );
}

// ************************************************************************* //
