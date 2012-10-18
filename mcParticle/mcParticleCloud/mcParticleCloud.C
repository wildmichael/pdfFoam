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

#include "mcParticleCloud.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fixedValueFvPatchField.H"
#include "boundBox.H"
#include "fvc.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "compressible/LES/LESModel/LESModel.H"
#include "mcProcessorBoundary.H"
#include "gradInterpolationConstantTet.H"
#include "uniqueOrder_FIX.H"

// * * * * * * * * * * * * * Local Helper Functions  * * * * * * * * * * * * //

namespace // anonymous
{

const Foam::compressible::turbulenceModel* getTurbulenceModel
(
    const Foam::objectRegistry& obr
)
{
    if (obr.foundObject<Foam::compressible::turbulenceModel>("RASProperties"))
    {
        return &obr.lookupObject<Foam::compressible::RASModel>
        (
            "RASProperties"
        );
    }
    else if (obr.foundObject<Foam::compressible::LESModel>("LESProperties"))
    {
        return &obr.lookupObject<Foam::compressible::LESModel>
        (
            "LESProperties"
        );
    }
    else
    {
        Foam::FatalErrorIn("getTurbulenceModel(const Foam::objectRegistry&)")
            << "No valid model for TKE calculation."
            << Foam::exit(Foam::FatalError);
        return 0;
    }
}


void constrainParticle
(
    Foam::mcParticleCloud& cloud,
    Foam::scalar dt,
    Foam::mcParticle& p
)
{
    using namespace Foam;
    const polyMesh& mesh = cloud.mesh();
    vector& u = p.Utracking();
    meshTools::constrainDirection(mesh, mesh.solutionD(), u);

    if (cloud.isAxiSymmetric())
    {
        point destPos = p.position() + dt*u;
        vector n = cloud.axis() ^ destPos;
        n /= mag(n);
        tensor T = rotationTensor(n, cloud.centrePlaneNormal());
        p.transformProperties(T);
        destPos = transform(T, destPos);
        // constrain to kill numerical artifacts
        meshTools::constrainDirection(mesh, mesh.geometricD(), destPos);
        // constrained tracking velocity to destPos
        u = (destPos - p.position())/dt;
    }
}


//- @todo This is a hack to work around annoying bug in OpenFOAM < 2.0
template<class DF>
void readIfPresent(DF& df)
{
    using namespace Foam;
    if
    (
        (df.readOpt() == IOobject::READ_IF_PRESENT && df.headerOk())
     || df.readOpt() == IOobject::MUST_READ
    )
    {
        df.readField(dictionary(df.readStream(df.typeName)), "value");
    }
}

const Foam::dimensionedScalar SMALL_MASS
(
    "SMALL_MASS",
    Foam::dimensionSet(1.,0.,0.,0.,0.,0.,0.),
    Foam::SMALL
);


const Foam::dimensionedScalar SMALL_VOLUME
(
    "SMALL_VOLUME",
    Foam::dimensionSet(0.,3.,0.,0.,0.,0.,0.),
    Foam::SMALL
);

} // anonymous namespace

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineParticleTypeNameAndDebug(mcParticle, 0);
    defineTemplateTypeNameAndDebug(Cloud<mcParticle>, 0);

}  // namespace Foam


inline bool Foam::mcParticleCloud::less::operator()
(
    mcParticle* one,
    mcParticle* two
) const
{
    return (one->eta()*one->m()) < (two->eta()*two->m());
}


inline bool Foam::mcParticleCloud::more::operator()
(
    mcParticle* one,
    mcParticle* two
) const
{
    return (one->eta()*one->m()) > (two->eta()*two->m());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcParticleCloud::mcParticleCloud
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cloudName,
    const compressible::turbulenceModel* turbModel,
    const volVectorField* U,
    const volScalarField* p,
    volScalarField* rho
)
:
    Cloud<mcParticle>(mesh, cloudName, false),
    mesh_(mesh),
    thermoDict_(dict),
    solutionDict_(mesh),
    runTime_(mesh.time()),
    deltaT_("particleDeltaT", dimTime, 100*SMALL),
    turbModel_
    (
        turbModel ? *turbModel : *getTurbulenceModel(mesh)
    ),
    Ufv_
    (
        U ? *U : mesh_.lookupObject<volVectorField>
                     (thermoDict_.lookupOrDefault<word>("UName", "U"))
    ),
    pfv_
    (
        p ? *p : mesh_.lookupObject<volScalarField>
                     (thermoDict_.lookupOrDefault<word>("pName", "p"))
    ),
    random_(55555+12345*Pstream::myProcNo()),
    scalarNames_(0),
    Nc_(mesh_.nCells()),
    histNp_(size()),

    cellParticleAddr_(Nc_),

    PaNIC_
    (
        IOobject
        (
            "PaNIC",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("PaNIC", dimless, 0)
    ),

    mMom_
    (
        IOobject
        (
            "mMoment",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mMoment", dimMass, 0)
    ),

    VMom_
    (
        IOobject
        (
            "VMoment",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("VMoment", dimVolume, 0)
    ),

    UMom_
    (
        IOobject
        (
            "UMoment",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("UMoment", dimMass*dimVelocity, vector::zero)
    ),

    PhiMom_(0),

    PhiPhiMom_(0),

    UUMom_
    (
        IOobject
        (
            "UUMoment",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("UUMoment", dimEnergy, symmTensor::zero)
    ),

    ownedScalarFields_(),

    rhocPdf_
    (
        rho ? *rho : const_cast<volScalarField&>(
            mesh_.lookupObject<volScalarField>(
                thermoDict_.lookupOrDefault<word>("rhoName", "rho")))
    ),

    rhocPdfInst_
    (
        IOobject
        (
            "rhoCloudPDFInst",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimDensity,
        rhocPdf_,
        rhocPdf_.boundaryField()
    ),

    pndcPdfInst_
    (
        IOobject
        (
            "pndCloudPDFInst",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimDensity,
        // Note: for axi-symmetric cases to *V/volumeOrArea later
        mMom_/mesh_.V(),
        rhocPdf_.boundaryField()
    ),

    pndcPdf_
    (
        IOobject
        (
            "pndCloudPDF",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimDensity,
        // Note: for axi-symmetric cases to *V/volumeOrArea later
        mMom_/mesh_.V(),
        rhocPdf_.boundaryField()
    ),

    UcPdf_
    (
        IOobject
        (
            "UCloudPDF",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimVelocity,
        UMom_/max(mMom_, SMALL_MASS),
        // Use the boundary conditions of U (FV)
        Ufv_.boundaryField()
    ),

    TaucPdf_
    (
        IOobject
        (
            "TauCloudPDF",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        turbModel_.R()
    ),

    kcPdf_
    (
        IOobject
        (
            "kCloudPDF",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimVelocity*dimVelocity,
        0.5*tr(TaucPdf_.dimensionedInternalField()),
        // Use the boundary conditions for k (FV)
        kfv()().boundaryField()
    ),

    PhicPdf_(),
    velocityModel_(),
    OmegaModel_(),
    mixingModel_(),
    reactionModel_(),
    positionCorrection_(),
    isAxiSymmetric_(false),
    CourantCoeffs_
    (
        IOobject
        (
            "CourantCoeffs_",
            mesh_.pointsInstance(),
            mesh_
        ),
        mesh_.surfaceInterpolation::deltaCoeffs()*mesh_.Sf()/mesh_.magSf()
    ),
    lostParticles_(*this),
    lostMass_(mesh_.V().size()),
    hNum_(0),
    deltaMass_
    (
        IOobject
        (
            "deltaMass",
            runTime_.timeName(),
            "uniform"/cloud::prefix/name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    massIn_
    (
        IOobject
        (
            "massIn",
            runTime_.timeName(),
            "uniform"/cloud::prefix/name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    massOut_
    (
        IOobject
        (
            "massOut",
            runTime_.timeName(),
            "uniform"/cloud::prefix/name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    cumDeltaMass_
    (
        IOobject
        (
            "cumDeltaMass",
            runTime_.timeName(),
            "uniform"/cloud::prefix/name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    cumMassIn_
    (
        IOobject
        (
            "cumMassIn",
            runTime_.timeName(),
            "uniform"/cloud::prefix/name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    cumMassOut_
    (
        IOobject
        (
            "cumMassOut",
            runTime_.timeName(),
            "uniform"/cloud::prefix/name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        0
    )
{
    bound(kcPdf_, solutionDict_.kMin());

    // HACK work around annoying bug in OpenFOAM < 2.0
    readIfPresent(mMom_);
    readIfPresent(VMom_);
    readIfPresent(UMom_);
    readIfPresent(UUMom_);

    CourantCoeffs_.boundaryField() /= 2.;

    initScalarFields();

    // Now that the fields exist, create the persisten models
    velocityModel_ = mcVelocityModel::New(*this, mesh_);
    OmegaModel_ = mcOmegaModel::New(*this, mesh_);
    mixingModel_ = mcMixingModel::New(*this, mesh_);
    reactionModel_ = mcReactionModel::New(*this, mesh_);
    positionCorrection_ = mcPositionCorrection::New(*this, mesh_);
    localTimeStepping_ = mcLocalTimeStepping::New(*this, mesh_);

    // Now determine whether this is an axi-symmetric case
    label nAxiSymmetric = 0;
    if (mesh_.nGeometricD() <= 2)
    {
        forAll(mesh_.boundaryMesh(), patchi)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchi];
            if (isA<wedgePolyPatch>(patch))
            {
                if (!isAxiSymmetric_)
                {
                    isAxiSymmetric_ = true;
                    const wedgePolyPatch& wpp = static_cast<const wedgePolyPatch&>(patch);
                    axis_ = wpp.axis();
                    centrePlaneNormal_ = wpp.centreNormal();
                    openingAngle_ =
                        2.*acos(wpp.patchNormal()&centrePlaneNormal_);
                    area_.reset(new DimensionedField<scalar, volMesh>
                        (
                            IOobject
                            (
                                "mcParticleCloud::area_",
                                runTime_.constant(),
                                mesh_,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                            ),
                            mesh_,
                            dimensionedScalar("A", dimVolume, 0)
                        ));
                    forAll(wpp, faceI)
                    {
                        area_()[wpp.faceCells()[faceI]] =
                            mag(wpp.faceAreas()[faceI]);
                    }
                }
                ++nAxiSymmetric;
            }
        }
        if (isAxiSymmetric_ && nAxiSymmetric != 2)
        {
            FatalErrorIn
            (
                "mcParticleCloud::mcParticleCloud"
                "("
                "    const fvMesh&,"
                "    const dictionary&,"
                "    const word&,"
                "    const compressible::turbulenceModel*,"
                "    const volVectorField*,"
                "    volScalarField*"
                ")"
            )  << "Only one pair of wedge patches allowed" << endl
               << exit(FatalError);
        }
    }

    initHNum();
    initBCHandlers();

    // Take care of statistical moments (make sure they are consistent)
    checkMoments();

    // Populate cloud
    if (returnReduce(size() > 0, andOp<bool>())) // if particle data found
    {
        mcParticle::readFields(*this);
    }
    else
    {
        Info<< "I am releasing particles initially!" << endl;
        initReleaseParticles();
    }
}


// If moments are not read correctly, initialize them.
void Foam::mcParticleCloud::checkMoments()
{
    bool readOk =
        PaNIC_.headerOk() &&
        mMom_.headerOk() &&
        VMom_.headerOk() &&
        UMom_.headerOk() &&
        UUMom_.headerOk() &&
        deltaMass_.headerOk() &&
        massIn_.headerOk() &&
        massOut_.headerOk() &&
        cumMassIn_.headerOk() &&
        cumMassOut_.headerOk();
    // Create moment fields
    label nPhi = PhicPdf_.size();
    PhiMom_.setSize(nPhi);
    PhiPhiMom_.setSize(PhiPhicPdf_.size());
    label PhiPhiI = 0;
    forAll(PhicPdf_, PhiI)
    {
        // Figure out a field name
        word name = PhicPdf_[PhiI]->name() + "Moment";
        dimensionSet dims = dimMass*PhicPdf_[PhiI]->dimensions();
        PhiMom_.set(PhiI, new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    name,
                    runTime_.timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(name, dims, 0)
            ));
        readIfPresent(PhiMom_[PhiI]);
        if (!PhiMom_[PhiI].headerOk())
        {
            readOk = false;
        }
        for (label PhiJ = PhiI; PhiJ != nPhi; ++PhiJ, ++PhiPhiI)
        {
            word PhiPhiName = PhiPhicPdf_[PhiPhiI]->name() + "Moment";
            dimensionSet dims =
                dimMass*PhicPdf_[PhiI]->dimensions()
               *PhicPdf_[PhiJ]->dimensions();
            PhiPhiMom_.set(PhiPhiI, new DimensionedField<scalar, volMesh>
                (
                    IOobject
                    (
                        PhiPhiName,
                        runTime_.timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar(PhiPhiName, dims, 0)
                ));
            readIfPresent(PhiPhiMom_[PhiPhiI]);
            if (!PhiPhiMom_[PhiPhiI].headerOk())
            {
                readOk = false;
            }
        }
    }
    if (returnReduce(readOk, andOp<bool>()))
    {
        Info<< "Moments read correctly." << endl;
        if (isAxiSymmetric_)
        {
            pndcPdf_.internalField() *= mesh_.V()/volumeOrArea();
            pndcPdfInst_.internalField() *= mesh_.V()/volumeOrArea();
        }
    }
    else
    {
        Info<< "Moments are missing. Forced re-initialization." << endl;
        initMoments();
    }
}


// Perform particle averaging to obtained cell-based values.
void Foam::mcParticleCloud::updateCloudPDF(scalar existWt)
{
    DimensionedField<scalar, volMesh> mMomInstant
    (
        IOobject
        (
            "mMomInstant",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("mMomInstant", dimMass, 0.0)
    );
    DimensionedField<scalar, volMesh> VMomInstant
    (
        IOobject
        (
            "VMomInstant",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("VMomInstant", dimVolume, 0.0)
    );
    DimensionedField<vector, volMesh> UMomInstant
    (
        IOobject
        (
            "UMomInstant",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("UMomInstant", dimMass*dimVelocity, vector::zero)
    );
    PtrList<DimensionedField<scalar, volMesh> >
        PhiMomInstant(PhiMom_.size()), PhiPhiMomInstant(PhiPhiMom_.size());
    label PhiPhiI = 0;
    forAll(PhiMom_, PhiI)
    {
        word name = PhiMom_[PhiI].name()+"Instant";
        PhiMomInstant.set
        (
            PhiI,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    name,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(name, PhiMom_[PhiI].dimensions(), 0.0)
            )
        );
        for (label PhiJ = PhiI; PhiJ != PhiMom_.size(); ++PhiJ, ++PhiPhiI)
        {
            word name = PhiPhiMom_[PhiPhiI].name()+"Instant";
            PhiPhiMomInstant.set
            (
                PhiPhiI,
                new DimensionedField<scalar, volMesh>
                (
                    IOobject
                    (
                        name,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar
                    (
                        name,
                        PhiPhiMom_[PhiPhiI].dimensions(),
                        0.0
                    )
                )
            );
        }
    }
    DimensionedField<symmTensor, volMesh> UUMomInstant
    (
        IOobject
        (
            "UUMomInstant",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("UUMomInstant", dimEnergy, symmTensor::zero)
    );

    PaNIC_ = 0;

    // Loop through particles to accumulate moments (0, 1, 2 order)
    // as well as particle number
    forAllConstIter(mcParticleCloud, *this, pIter)
    {
        const mcParticle& p = pIter();
        label cellI = p.cell();
        const vector& Up = p.UParticle();

        ++PaNIC_[cellI];

        const scalar mpd = p.eta()*massPerDepth(p);
        mMomInstant[cellI] += mpd;
        VMomInstant[cellI] += mpd/p.rho();
        UMomInstant[cellI] += mpd*p.UParticle();
        PhiPhiI = 0;
        forAll(PhicPdf_, PhiI)
        {
            PhiMomInstant[PhiI][cellI] += mpd*p.Phi()[PhiI];
            for (label PhiJ = PhiI; PhiJ != PhicPdf_.size(); ++PhiJ, ++PhiPhiI)
            {
                PhiPhiMomInstant[PhiPhiI][cellI] +=
                    mpd*p.Phi()[PhiI]*p.Phi()[PhiJ];
            }
        }
        UUMomInstant[cellI] += mpd*symm(Up*Up);
    }

    scalar newWt = 1.0 - existWt;
    // Do time-averaging of moments and compute mean fields
    mMom_  = existWt * mMom_  + newWt * mMomInstant;
    DimensionedField<scalar, volMesh> mMomBounded = max(mMom_, SMALL_MASS);
    pndcPdfInst_.internalField() = mMomInstant/volumeOrArea();
    pndcPdfInst_.correctBoundaryConditions();
    pndcPdf_.internalField() = mMom_/volumeOrArea();
    pndcPdf_.correctBoundaryConditions();

    VMom_  = existWt * VMom_  + newWt * VMomInstant;
    DimensionedField<scalar, volMesh> VMomBounded = max(VMom_, SMALL_VOLUME);
    rhocPdfInst_.internalField() = mMomInstant/max(VMomInstant, SMALL_VOLUME);
    rhocPdfInst_.correctBoundaryConditions();
    rhocPdf_.internalField()   = mMom_ / VMomBounded;
    rhocPdf_.correctBoundaryConditions();

    UMom_  = existWt * UMom_  + newWt * UMomInstant;
    UcPdf_.internalField()   = UMom_ / mMomBounded;
    UcPdf_.correctBoundaryConditions();

    PhiPhiI = 0;
    forAll(PhicPdf_, PhiI)
    {
        PhiMom_[PhiI] = existWt * PhiMom_[PhiI] + newWt * PhiMomInstant[PhiI];
        PhicPdf_[PhiI]->internalField() = PhiMom_[PhiI] / mMomBounded;
        PhicPdf_[PhiI]->correctBoundaryConditions();
        for (label PhiJ = PhiI; PhiJ != PhicPdf_.size(); ++PhiJ, ++PhiPhiI)
        {
            PhiPhiMom_[PhiPhiI] =
                existWt*PhiPhiMom_[PhiPhiI] + newWt*PhiPhiMomInstant[PhiPhiI];
            PhiPhicPdf_[PhiPhiI]->internalField() =
                PhiPhiMom_[PhiPhiI]/mMomBounded
              - (
                    PhicPdf_[PhiI]->dimensionedInternalField()
                   *PhicPdf_[PhiJ]->dimensionedInternalField()
                );
            PhiPhicPdf_[PhiPhiI]->correctBoundaryConditions();
        }
    }

    UUMom_ = existWt*UUMom_ + newWt*UUMomInstant;
    TaucPdf_.internalField() =
        (
            UUMom_/mMomBounded
          - symm(UcPdf_*UcPdf_)().dimensionedInternalField()
        );
    TaucPdf_.correctBoundaryConditions();

    kcPdf_.internalField()   = 0.5 * tr(TaucPdf_.internalField());
    kcPdf_.correctBoundaryConditions();
    bound(kcPdf_, solutionDict_.kMin());
}


// Update particle-cell addressing: for each cell, what are the addresses of
// the particles I host?
void Foam::mcParticleCloud::particleNumberControl()
{
    if (!solutionDict_.enableParticleNumberControl()) return;

    List<cellPopStatus> cellPopFlag(Nc_, NORMAL);
    const label Npc = solutionDict_.particlesPerCell();
    const scalar cloneAt = solutionDict_.cloneAt();
    const scalar eliminateAt = solutionDict_.eliminateAt();

    forAll(PaNIC_, celli)
    {
        label np = round(PaNIC_[celli]);

        // classify the particle # health condition of each cell
        if (np < 1)
        {
            cellPopFlag[celli] = EMPTY;
        }
        else if (np < round(Npc*cloneAt))
        {
            cellPopFlag[celli] = TOOFEW;
        }
        else if (np > round(Npc*eliminateAt))
        {
            cellPopFlag[celli] = TOOMANY;
        }

        // clear old list
        cellParticleAddr_[celli].clear();
    }

    labelList ncpi(Nc_, 0);

    forAllIter(mcParticleCloud, *this, pIter)
    {
        mcParticle & p = pIter();
        label celli = p.cell();
        // If a mcParticleList is necessary for this cell, construct it.
        if (cellPopFlag[celli] > NORMAL) // too few or too many particles
        {
            // according to ascending order of mass (if too many particles)
            // or descending order of mass (if too few particle)
            mcParticleList & cepl = cellParticleAddr_[celli];
            if (cepl.size() < 1)
            {
                cepl.setSize(round(PaNIC_[celli]), NULL);
            }
            cepl[ncpi[celli]++] = &p;
        }
    }

    // Sort the lists with particles, and perform particle number control
    forAll(cellPopFlag, celli)
    {
        if (cellPopFlag[celli] == TOOFEW)
        {
            cloneParticles(celli);
        }
        else if ( cellPopFlag[celli] == TOOMANY )
        {
            eliminateParticles(celli);
        }
    }

#if 0
    // Debug only: check the list
    if (debug)
    {
        forAll(cellPopFlag, celli)
        {
            if (cellPopFlag[celli] == TOOFEW || cellPopFlag[celli] == TOOMANY)
            {
                Info<< "size is " << cellParticleAddr_[celli].size()
                    << ", flag is " << cellPopFlag[celli]
                    << " := " << endl;
                mcParticleList & cepl = cellParticleAddr_[celli];

                forAll(cepl, pci)
                {
                    if (cepl[pci])
                    {
                        Info<< "ID= " << cepl[pci]->origId() << ", "
                            << "m= " << cepl[pci]->m() << endl;
                    }
                }
            }
        }
    }
#endif
}


// Split the n heaviest particles
void Foam::mcParticleCloud::cloneParticles(label celli)
{
    // no. particle to reproduce
    label n = solutionDict_.particlesPerCell() - round(PaNIC_[celli]);
    n = min(round(PaNIC_[celli]), n);

    sort(cellParticleAddr_[celli], more());

    vectorList positions = randomPointsInCell(n, celli);

    for (label particleI=0; particleI < n; particleI++)
    {
        mcParticle& p = *(cellParticleAddr_[celli][particleI]);
        // Half my mass
        p.m() /= 2.0;
        // create a new particle like myself
        autoPtr<mcParticle> ptrNew = p.clone();
        ptrNew().position() = positions[particleI];

        addParticle(ptrNew.ptr());
    }
    PaNIC_[celli] += n;

    histNp_ += n;
    if (debug)
    {
        Pout<< "Cloned " << n << " particles in cell " << celli << endl;
    }
}


// As name suggests
void Foam::mcParticleCloud::eliminateParticles(label celli)
{
    label ncur = round(PaNIC_[celli]);
    // no. particle to eliminate
    label nx =  ncur - solutionDict_.particlesPerCell();

    // All particles in celli
    mcParticleList& popAll = cellParticleAddr_[celli];
    // The two randomly selected populations
    SLList<mcParticle*> popA, popB;
    // Masses of the random populations
    scalar mA = 0., mB = 0.;

    // Probability to select a particle into any of the random populations
    scalar P = (2.*nx)/ncur;

    // Create two populations of size nx (in average)
    forAllIter(mcParticleList, popAll, pIter)
    {
        if (random().scalar01() < P)
        {
            mcParticle* p = *pIter;
            scalar meta = p->eta()*p->m();
            if (random().scalar01() < 0.5)
            {
                popA.append(p);
                mA += meta;
            }
            else
            {
                popB.append(p);
                mB += meta;
            }
        }
    }

    // Only continue if one of the populations has some mass
    if (mA > 0 || mB > 0)
    {
        // Randomly pick one of the populations for deletion
        // (proportional to relative mass of the other population)
        SLList<mcParticle*> *popDel, *popKeep;
        P = mB/(mA + mB);
        scalar s;
        if (random().scalar01() < P)
        {
            s = 1./P;
            popDel = &popA;
            popKeep = &popB;
        }
        else
        {
            s = (mA + mB)/mA;
            popDel = &popB;
            popKeep = &popA;
        }

        // Scale masses of particles in popKeep
        forAllIter(SLList<mcParticle*>, *popKeep, pIter)
        {
            (**pIter).m() *= s;
        }

        // Delete particles in popDel
        label nKilled = 0;
        forAllIter(SLList<mcParticle*>, *popDel, pIter)
        {
            ++nKilled;
            deleteParticle(**pIter);
        }
        PaNIC_[celli] -= nKilled;

        if (debug)
        {
            Pout<< "Eliminated " << nKilled << " of " << ncur
                << " particles in cell " << celli << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::mcParticleCloud::evolve()
{

    // Correct boundary conditions
    forAll(boundaryHandlers_, boundaryI)
    {
        mcBoundary& b = boundaryHandlers_[boundaryI];
        b.massIn() = 0.;
        b.massOut() = 0.;
        b.correct(false);
    }
    // Integrate scalars across domain and reset correction velocity
    scalarField deltaMassInst(deltaMass_.size(), 0.);
    forAllIter(mcParticleCloud, *this, pIter)
    {
        mcParticle& p = pIter();
        scalar mpd = massPerDepth(p);
        deltaMassInst[0] -= mpd;
        forAll(conservedScalars_, csI)
        {
            deltaMassInst[csI+1] -= mpd*p.Phi()[conservedScalars_[csI]];
        }
        p.Ucorrection() = vector::zero;
    }
    positionCorrection_().correct();
    // Reset lost mass
    lostMass_ = 0;

    // First half-step
    //////////////////

    forAllIter(mcParticleCloud, *this, pIter)
    {
        mcParticle& p = pIter();
        p.nSteps() = 0;
        p.Utracking() = p.UParticle() + p.Ucorrection();
        constrainParticle(*this, deltaT_.value()/2, p);
        p.reflected() = false;
        // Store old data
        p.UParticleOld() = p.UParticle();
        p.positionOld() = p.position();
        p.cellOld() = p.cell();
        p.faceOld() = p.face();
        p.procOld() = Pstream::myProcNo();
    }

    mcParticle::trackData td1(*this, deltaT_.value()/2.);
    Cloud<mcParticle>::move(td1);

    // Evaluate models at deltaT/2
    forAllIter(mcParticleCloud, *this, pIter)
    {
        pIter().nSteps() = 0;
        computeCourantNo(pIter());
    }
    OmegaModel_().correct();
    mixingModel_().correct();
    reactionModel_().correct();
    velocityModel_().correct();
    localTimeStepping_().correct();

    // Send back to original processor
    if (Pstream::parRun())
    {
        List<IDLList<mcParticle> > transferList(Pstream::nProcs());
        forAllIter(mcParticleCloud, *this, pIter)
        {
            mcParticle& p = pIter();

            // If this is a parallel run and the particle switch processor in
            // the first half step, put it in the transfer list
            if (Pstream::parRun() && p.procOld() != Pstream::myProcNo())
            {
                transferList[p.procOld()].append(this->remove(&p));
            }
        }
        // List of the numbers of particles to be transfered across the
        // processor patches
        labelList nsTransPs(transferList.size());

        forAll(transferList, i)
        {
            nsTransPs[i] = transferList[i].size();
        }

        // List of the numbers of particles to be transfered across the
        // processor patches for all the processors
        labelListList allNTrans(Pstream::nProcs());
        allNTrans[Pstream::myProcNo()] = nsTransPs;
        combineReduce(allNTrans, combineNsTransPs());

        forAll(transferList, i)
        {
            if (transferList[i].size())
            {
                OPstream particleStream(Pstream::blocking, i);
                particleStream << transferList[i];
            }
        }

        forAll(allNTrans, i)
        {
            label nRecPs = allNTrans[i][Pstream::myProcNo()];

            if (nRecPs)
            {
                IPstream particleStream(Pstream::blocking, i);
                IDLList<mcParticle> newParticles
                (
                    particleStream,
                    mcParticle::iNew(*this)
                );

                forAllIter(IDLList<mcParticle>, newParticles, newpIter)
                {
                    mcParticle& newp = newpIter();
                    addParticle(newParticles.remove(&newp));
                }
            }
        }
    }

    // Estimate particle velocity as 0.5*(U^{n}+U^{n+1}) and put particles back
    // to their original position. For particles that have been reflected,
    // decay to first-order integration.
    forAllIter(mcParticleCloud, *this, pIter)
    {
        mcParticle& p = pIter();

        p.position() = p.positionOld();
        p.cell() = p.cellOld();
        p.face() = p.faceOld();

        if (p.reflected())
        {
            p.Utracking() = p.UParticleOld();
        }
        else
        {
            p.Utracking() = 0.5*(p.UParticleOld() + p.UParticle());
        }

        // Add numerical diffusion (random walk)
        const scalar& DNum = solutionDict_.DNum().value();
        scalar C =
            DNum*sqrt(hNum_[p.cell()]*deltaT_.value()*mag(p.Utracking()))
           /deltaT_.value();
        p.Utracking() +=
            vector
            (
                C*random_.GaussNormal(),
                C*random_.GaussNormal(),
                C*random_.GaussNormal()
            );

        // Add correction velocity
        p.Utracking() += p.Ucorrection();

        constrainParticle(*this, deltaT_.value(), p);
    }

    // Second half-step
    //////////////////

    mcParticle::trackData td2(*this, deltaT_.value());
    Cloud<mcParticle>::move(td2);

    // Correct boundary conditions
    scalarField massInInst(massIn_.size(), 0.);
    scalarField massOutInst(massOut_.size(), 0.);
    forAll(boundaryHandlers_, boundaryI)
    {
        mcBoundary& b = boundaryHandlers_[boundaryI];
        b.correct(true);
        massInInst += b.massIn();
        massOutInst += b.massOut();
    }

    // Extract statistical averaging to obtain mesh-based quantities
    const dimensionedScalar& avgTime = solutionDict_.averagingTime();
    scalar existWt = 1.0/(1.0 + (deltaT_/(avgTime - deltaT_)).value());
    updateCloudPDF(existWt);

    particleNumberControl();

    // Redistribute lost mass and reset lost mass counter
    // FIXME Doing this after the extraction is probably suboptimal
    lostMass_ = lostMass_ / max(PaNIC_.internalField(), 1.);
    forAllIter(mcParticleCloud, *this, pIter)
    {
        mcParticle& p = pIter();
        p.m() += lostMass_[p.cell()];
    }
    lostMass_ = 0;

    if (debug)
    {
        assertPopulationHealth();
    }

    scalarField rhoErr =
        (
            pndcPdf_.internalField()
          - rhocPdf_.internalField()
        )/rhocPdf_.internalField();
    scalar rhoRes = gAverage(mag(rhoErr));
    Info<< "Cloud " << name() << nl
        << "    deltaT = " << deltaT_.value() << nl
        << "    pmd relative error: averageMag = " << rhoRes
        << ", min = " << gMin(rhoErr)
        << ", max = " << gMax(rhoErr)
        << nl;

    label nLostPart = returnReduce(lostParticles_.size(), sumOp<label>());
    Info<< "    number of lost particles: " << nLostPart << nl;
#ifdef FULLDEBUG
    volVectorField stabUfv = Ufv_;
    forAll(stabUfv, cellI)
    {
        stabUfv[cellI] = stabilise(stabUfv[cellI], 100*SMALL);
    }
    forAll(stabUfv.boundaryField(), patchI)
    {
        forAll(stabUfv.boundaryField()[patchI], faceI)
        {
            vector& Up = stabUfv.boundaryField()[patchI][faceI];
            Up = stabilise(Up, 100*SMALL);
        }
    }
    volScalarField deltaU
    (
        IOobject
        (
            "deltaU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag
        (
            cmptDivide
            (
                UcPdf_ - Ufv_,
                stabUfv
            )
        )
    );
    if (mesh_.time().outputTime())
    {
        deltaU.write();
    }
    Info<< "DEBUG: maximum deltaU = " << gMax(deltaU) << endl;
#endif
    // Mass and integrated scalars after the evolution done, compute maximum
    // Courant number
    scalar totalMass = fvc::domainIntegrate(rhocPdfInst_).value();
    scalar totalParticleMass = fvc::domainIntegrate(pndcPdfInst_).value();
    scalar UMax = 0, UcorrMax = 0, etaMax = 0, etaMin = HUGE;
    scalar CoMax = 0;
    forAllConstIter(mcParticleCloud, *this, pIter)
    {
        const mcParticle& p = pIter();
        scalar mpd = massPerDepth(p);
        deltaMassInst[0] += mpd;
        forAll(conservedScalars_, csI)
        {
            deltaMassInst[csI+1] += mpd*p.Phi()[conservedScalars_[csI]];
        }
        UMax = max(UMax, mag(p.UParticle()));
        UcorrMax = max(UcorrMax, mag(p.Ucorrection()));
        etaMax = max(etaMax, p.eta());
        etaMin = min(etaMin, p.eta());
        CoMax = max(CoMax, p.Co());
    }
    reduce(UMax, maxOp<scalar>());
    reduce(UcorrMax, maxOp<scalar>());
    reduce(etaMin, maxOp<scalar>());
    reduce(etaMax, maxOp<scalar>());
    reduce(deltaMassInst,  sumOp<scalarField>());
    reduce(massInInst,  sumOp<scalarField>());
    reduce(massOutInst, sumOp<scalarField>());
    reduce(CoMax, maxOp<scalar>());

    static bool firstRun = true;
    if (!firstRun)
    {
        deltaMass_  = 0.999*deltaMass_  + 0.001*deltaMassInst;
        massIn_  = 0.999*massIn_  + 0.001*massInInst;
        massOut_ = 0.999*massOut_ + 0.001*massOutInst;

        cumDeltaMass_ += deltaMassInst;
        cumMassIn_ += massInInst;
        cumMassOut_ += massOutInst;

        scalarList massErr = (deltaMass_ + massIn_ + massOut_)/stabilise(massIn_, SMALL);
        scalarList massErrInst = (deltaMassInst + massInInst + massOutInst)/stabilise(massInInst, SMALL);
        scalarList cumMassErr = (cumDeltaMass_ + cumMassIn_ + cumMassOut_)/stabilise(cumMassIn_, SMALL);

        Info<< "    instant mass:  density = " << totalMass
            << ", particle = " << totalParticleMass
            << ", error/densityMass = "
            << (totalParticleMass - totalMass)/totalMass << nl
            << "    min(rhoCloudPdfInst) = " << gMin(rhocPdfInst_)
            << ", max(rhoCloudPdfInst) = " << gMax(rhocPdfInst_) << nl
            << "    min(pndCloudPdfInst) = " << gMin(pndcPdfInst_)
            << ", max(pndCloudPdfInst) = " << gMax(pndcPdfInst_) << nl;

        Info<< "    particle mass flux error: "
            << "instant = " << massErrInst[0]
            << ", mean = " << massErr[0]
            << ", cumulative = " << cumMassErr[0] << nl;

        forAll(conservedScalars_, csI)
        {
            label i = conservedScalars_[csI];
            Info<< "    " << scalarNames_[i] << " flux error: "
                << "instant = " << massErrInst[csI+1]
                << ", mean = " << massErr[csI+1]
                << ", cumulative = " << cumMassErr[csI+1] << nl;
        }
    }
    firstRun = false;

    Info<< "    max(UParticle) = " << UMax
        << ", max(Ucorrection) = " << UcorrMax
        << ", min(eta) = " << etaMin
        << ", max(eta) = " << etaMax
        << nl;

    // Finally, update deltaT for next time step
    deltaT_.value() = solutionDict_.CFL()/CoMax;

    return rhoRes;
}


// Initialization: populate the FV field with particles
void Foam::mcParticleCloud::initReleaseParticles()
{
    clear();
    const label Npc = solutionDict_.particlesPerCell();
    // Populate each cell with particlesPerCell() particles in each cell
    forAll(Ufv_, celli)
    {
        scalar m = mesh_.V()[celli]*rhocPdf_[celli]/Npc;
        // TODO fv or pdf?
        vector Updf = Ufv_[celli];
        scalar urms = sqrt(2./3.*kfv()()[celli]);
        vector uscales(urms, urms, urms);
        scalarField Phi(PhicPdf_.size());
        forAll(Phi, PhiI)
        {
            Phi[PhiI] = (*PhicPdf_[PhiI])[celli];
        }

        particleGenInCell(celli, Npc, m, Updf, uscales, Phi);
    }
    // correct mass according to local time-stepping
    localTimeStepping_().correct();
    forAllIter(mcParticleCloud, *this, pIter)
    {
        pIter().m() /= pIter().eta();
    }
    OmegaModel_().correct();
    reactionModel_().correct();
    // writeFields();
}


// Randomly generate N particles in celli, with provided cell-based
// values and the scale of velocity fluctuation
void Foam::mcParticleCloud::particleGenInCell
(
    label celli,
    label N,
    scalar m,
    const vector& Updf,
    const vector& uscales,
    const scalarField& Phi,
    const vector& shift,
    label  ghost
)
{

    label Npgen = 0;
    DynamicList<mcParticle*> genParticles;
    genParticles.reserve(N);

    vectorList positions = randomPointsInCell(N, celli);

    for (label i=0; i<N; ++i)
    {
        // What is the distribution of u?
        // How to enforce the component-wise correlations < u_i, u_j >?
        // TODO generate correlated fluctuations
        vector u(
            random().GaussNormal() * uscales.x(),
            random().GaussNormal() * uscales.y(),
            random().GaussNormal() * uscales.z()
            );
        vector UParticle = u + Updf;

        mcParticle* ptr = new mcParticle
        (
            *this,
            positions[i],
            celli,
            m,
            UParticle,
            Phi,
            shift,
            ghost
        );

        addParticle(ptr);
        genParticles.append(ptr);
        ++Npgen;
    }

    adjustAxiSymmetricMass(genParticles);

    if (Npgen < N)
    {
        FatalErrorIn("mcParticleCloud::initReleaseParticles()")
            << "Only " << Npgen << " particles generated for cell "
            << celli << nl << "Something is wrong" << exit(FatalError);
    }
    else // update count
    {
        histNp_ += Npgen;
    }

}


// Initialize statistical moments from FV data
void Foam::mcParticleCloud::initMoments()
{
    mMom_  = rhocPdf_*volumeOrArea();
    rhocPdfInst_.internalField() = rhocPdf_.internalField();
    pndcPdf_.internalField() = rhocPdf_.internalField();
    pndcPdfInst_.internalField() = rhocPdf_.internalField();

    VMom_  = mMom_ / rhocPdf_;

    UMom_  = mMom_ * Ufv_;
    UcPdf_.internalField()   = Ufv_.internalField();
    UcPdf_.correctBoundaryConditions();

    label PhiPhiI = 0;
    forAll(PhicPdf_, PhiI)
    {
        PhiMom_[PhiI] = mMom_ * (*PhicPdf_[PhiI]);
        for (label PhiJ = PhiI; PhiJ != PhicPdf_.size(); ++PhiJ, ++PhiPhiI)
        {
            PhiPhiMom_[PhiPhiI] =
                mMom_
               *(
                   PhiPhicPdf_[PhiPhiI]->dimensionedInternalField()
                 + (
                       PhicPdf_[PhiI]->dimensionedInternalField()
                      *PhicPdf_[PhiJ]->dimensionedInternalField()
                   )
                );
        }
    }

    UUMom_ =
        mMom_
       *(
           turbulenceModel().R()()
         + symm(UcPdf_*UcPdf_)().dimensionedInternalField()
        );

    kcPdf_.internalField()   = turbulenceModel().k()().internalField();
    kcPdf_.correctBoundaryConditions();
    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>
        (
            thermoDict_.lookupOrDefault<word>("phiName", "phi")
        );
    deltaMass_.setSize(conservedScalars_.size()+1, 0.);
    massIn_.setSize(deltaMass_.size(), 0.);
    massOut_.setSize(deltaMass_.size(), 0.);
    forAll(phi.boundaryField(), patchI)
    {
        const scalarField& phiPatch = phi.boundaryField()[patchI];
        forAll(phiPatch, faceI)
        {
            scalarField& mass =
                phiPatch[faceI] < 0 ? massIn_ : massOut_;
            mass[0] += -phiPatch[faceI];
            forAll(conservedScalars_, csI)
            {
                const scalarField& PhiPatch =
                    PhicPdf_[conservedScalars_[csI]]->boundaryField()[patchI];
                mass[csI+1] += -PhiPatch[faceI]*phiPatch[faceI];
            }
        }
    }
    deltaMass_ = massIn_ + massOut_;
    cumDeltaMass_.setSize(deltaMass_.size(), 0.);
    cumMassIn_.setSize(deltaMass_.size(), 0.);
    cumMassOut_.setSize(deltaMass_.size(), 0.);
}


void Foam::mcParticleCloud::initScalarFields()
{
    // Find scalar fields
    if (thermoDict_.found("scalarFields"))
    {
        thermoDict_.lookup("scalarFields") >> scalarNames_;
        // Check for duplicates
        labelList order;
        uniqueOrder_FIX(scalarNames_, order);
        if (scalarNames_.size() != order.size())
        {
            FatalErrorIn
            (
                "mcParticleCloud::initScalarFields()"
            )
                << "The list "
                << thermoDict_.lookupEntry("scalarFields", false, false).name()
                << " contains duplicate entries.\n"
                << exit(FatalError);
        }
    }
    label nScalarFields = scalarNames_.size();
    PhicPdf_.setSize(nScalarFields);
    PhiPhicPdf_.setSize(label(nScalarFields*(nScalarFields + 1.)/2.));
    label PhiPhiI = 0;
    forAll(scalarNames_, fieldI)
    {
        if (mesh_.foundObject<volScalarField>(scalarNames_[fieldI]))
        {
            // Try to find that field
            PhicPdf_[fieldI] =  &const_cast<volScalarField&>
                (mesh_.lookupObject<volScalarField>(scalarNames_[fieldI]));
        }
        else
        {
            // Field doesn't exist already, so insert a new one into
            // ownedScalarFields_
            Info<< "Creating mcParticleCloud-owned field "
                << scalarNames_[fieldI] << endl;
            ownedScalarFields_.insert
            (
                new volScalarField
                (
                    IOobject
                    (
                        scalarNames_[fieldI],
                        runTime_.timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_
                )
            );
            PhicPdf_[fieldI] = &ownedScalarFields_.first();
        }
        // scalar covariance fields
        for (label fieldJ = fieldI; fieldJ != nScalarFields; ++fieldJ, ++PhiPhiI)
        {
            word name = scalarNames_[fieldI]+scalarNames_[fieldJ] + "Cov";
            if (mesh_.foundObject<volScalarField>(name))
            {
                // Try to find that field
                PhiPhicPdf_[PhiPhiI] =  &const_cast<volScalarField&>
                    (mesh_.lookupObject<volScalarField>(name));
            }
            else
            {
                // Field doesn't exist already, so insert a new one into
                // ownedScalarFields_
                Info<< "Creating mcParticleCloud-owned field "
                    << name << endl;
                ownedScalarFields_.insert
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            name,
                            runTime_.timeName(),
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh_
                    )
                );
                PhiPhicPdf_[PhiPhiI] = &ownedScalarFields_.first();
            }
        }
    }

    // Find labels of mixed scalars
    scalarNames_.setSize(nScalarFields);
    wordList mixedScalarNames;
    if (thermoDict_.found("mixedScalars"))
    {
        thermoDict_.lookup("mixedScalars") >> mixedScalarNames;
        // Check for duplicates
        labelList order;
        uniqueOrder_FIX(mixedScalarNames, order);
        if (mixedScalarNames.size() != order.size())
        {
            FatalErrorIn
            (
                "mcParticleCloud::initScalarFields()"
            )
                << "The list "
                << thermoDict_.lookupEntry("mixedScalars", false, false).name()
                << " contains duplicate entries.\n"
                << exit(FatalError);
        }
    }
    label nMixedScalarFields = mixedScalarNames.size();
    mixedScalars_.setSize(nMixedScalarFields);
    forAll(mixedScalarNames, mixedScalarI)
    {
        bool found = false;
        forAll(scalarNames_, fieldI)
        {
            if (mixedScalarNames[mixedScalarI] == scalarNames_[fieldI])
            {
                mixedScalars_[mixedScalarI] = fieldI;
                found = true;
                break;
            }
        }
        if (!found)
        {
            FatalErrorIn
            (
                "mcParticleCloud::initScalarFields()"
            )
                << "No such scalar field: " << mixedScalarNames[mixedScalarI] << "\n"
                << "Available field names are:\n" << scalarNames_ << "\n"
                << exit(FatalError);
        }
    }
    // Find conserved scalars
    wordList conservedScalarNames;
    if (thermoDict_.found("conservedScalars"))
    {
        thermoDict_.lookup("conservedScalars") >> conservedScalarNames;
        // Check for duplicates
        labelList order;
        uniqueOrder_FIX(conservedScalarNames, order);
        if (conservedScalarNames.size() != order.size())
        {
            FatalErrorIn
            (
                "mcParticleCloud::initScalarFields()"
            )
                << "The list "
                << thermoDict_.lookupEntry("conservedScalars", false, false).name()
                << " contains duplicate entries.\n"
                << exit(FatalError);
        }
    }
    label nConservedScalarFields = conservedScalarNames.size();
    conservedScalars_.setSize(nConservedScalarFields);
    forAll(conservedScalarNames, conservedScalarI)
    {
        bool found = false;
        forAll(scalarNames_, fieldI)
        {
            if (conservedScalarNames[conservedScalarI] == scalarNames_[fieldI])
            {
                conservedScalars_[conservedScalarI] = fieldI;
                found = true;
                break;
            }
        }
        if (!found)
        {
            FatalErrorIn
            (
                "mcParticleCloud::initScalarFields()"
            )
                << "No such scalar field: " << conservedScalarNames[conservedScalarI] << "\n"
                << "Available field names are:\n" << scalarNames_ << "\n"
                << exit(FatalError);
        }
    }
}


void Foam::mcParticleCloud::initBCHandlers()
{
    boundaryHandlers_.clear();
    boundaryHandlers_.setSize(mesh_.boundaryMesh().size());
    const dictionary& bd = thermoDict_.subDict("boundaryHandlers");
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        autoPtr<mcBoundary> bcHandler;
        if (isA<processorPolyPatch>(pp))
        {
            bcHandler.reset(new mcProcessorBoundary
            (
                *this,
                patchI
            ));
        }
        else
        {
            bcHandler = mcBoundary::New
            (
                *this,
                patchI,
                bd.subDict(pp.name())
            );
        }
        boundaryHandlers_.set
        (
            patchI,
            bcHandler
        );
    }
}


void Foam::mcParticleCloud::initHNum()
{
    const pointField& points = mesh_.points();
    const faceList& faces = mesh_.faces();
    const cellList& cells = mesh_.cells();
    hNum_.setSize(cells.size());
    forAll(hNum_, celli)
    {
        labelList cellVertices = cells[celli].labels(faces);
        vector bbmax = -GREAT*vector::one;
        vector bbmin = GREAT*vector::one;
        forAll (cellVertices, vertexi)
        {
            bbmax = max(bbmax, points[cellVertices[vertexi]]);
            bbmin = min(bbmin, points[cellVertices[vertexi]]);
        }
        vector dx
            (
                bbmax.x()-bbmin.x(),
                bbmax.y()-bbmin.y(),
                bbmax.z()-bbmin.z()
            );
        if (mesh_.nGeometricD() < 3)
        {
            dx.x() = mesh_.geometricD()[0] < 0 ? 1. : dx.x();
            dx.y() = mesh_.geometricD()[1] < 0 ? 1. : dx.y();
            dx.y() = mesh_.geometricD()[2] < 0 ? 1. : dx.z();
            hNum_[celli] = sqrt(dx.x()*dx.y()*dx.z());
        }
        else
        {
            hNum_[celli] = cbrt(dx.x()*dx.y()*dx.z());
        }
    }
}


void Foam::mcParticleCloud::info() const
{
    Info<< "Cloud: " << this->name() << nl
        << "    Current number of particles         = "
        << returnReduce(this->size(), sumOp<label>()) << nl;
//      << "    Totally number of particles existed = "
//      << returnReduce(histNp_,      sumOp<label>()) << nl;
}


void Foam::mcParticleCloud::assertPopulationHealth() const
{
    Info<< "Asserting particle population health" << endl;
    bool cloudHealthy = size() > 0;
    OStringStream reasons;
    if (!cloudHealthy)
    {
       reasons
           << "\tNo particles in the cloud" << nl;
    }

    labelList nPart(Ufv_.size(), 0);
    forAllConstIter(mcParticleCloud, *this, pIter)
    {
        const mcParticle& p = pIter();
        nPart[p.cell()] += 1;
        if (p.cell() != mesh_.findCell(p.position()))
        {
            cloudHealthy = false;
            reasons
                << "\tParticle " << p.origId() << " is not in cell "
                << p.cell() << nl;
        }
        if (p.ghost())
        {
            cloudHealthy = false;
            reasons
                << "\tParticle " << p.origId()
                << " is a ghost particle" << nl;
        }
        if (mag(p.shift()) > 100*SMALL)
        {
            cloudHealthy = false;
            reasons
                << "\tParticle " << p.origId()
                << " has a non-zero shift vector" << nl;
        }
    }
    forAll(nPart, cellI)
    {
        if (!nPart[cellI])
        {
            cloudHealthy = false;
            reasons << "\tCell " << cellI << " has no particles" << nl;
        }
    }
    if (!cloudHealthy)
    {
        FatalErrorIn("mcParticleCloud::assertPopulationHealth()")
            << "Particle health check not passed:" << nl << reasons.str()
            << exit(FatalError);
    }

}


Foam::vectorList
Foam::mcParticleCloud::randomPointsInCell(label n, label celli)
{
    autoPtr<vectorList> ppoints(new vectorList(n, point::zero));
    vectorList& points = ppoints();

    boundBox cellbb
    (
        pointField
        (
            mesh_.points(),
            mesh_.cellPoints()[celli]
        ),
        false
    );

    vector minb = cellbb.min();
    vector dimb = cellbb.max() - minb;

    for (label i=0; i<n;)
    {
        // Relative coordinate [0, 1] in this cell
        vector xi = random().vector01();
        // Random offset from min point
        scalar rx = min(max(10.0*SMALL, xi.x()), 1.0-10.0*SMALL);
        scalar ry = min(max(10.0*SMALL, xi.y()), 1.0-10.0*SMALL);
        scalar rz = min(max(10.0*SMALL, xi.z()), 1.0-10.0*SMALL);
        vector offsetRnd(rx*dimb.x(), ry*dimb.y(), rz*dimb.z());

        // Generate a particle position
        vector position = minb + offsetRnd;

        // If the case has reduced dimensionality, put the coordinate of the
        // reduced dimension onto the coordinate plane
        if (mesh_.nGeometricD() <= 2)
        {
            meshTools::constrainDirection(mesh_, mesh_.geometricD(), position);
        }

        if (mesh_.pointInCell(position, celli))
        {
            points[i++] = position;
        }
    }
    return ppoints;
}


void Foam::mcParticleCloud::notifyLostParticle(const Foam::mcParticle& p)
{
    lostMass_[p.cell()] += p.m();
    lostParticles_.add(p);
}


void Foam::mcParticleCloud::adjustAxiSymmetricMass
(
    UList<mcParticle*>& particles
)
{
    if (!isAxiSymmetric_)
    {
        return;
    }
    scalarList r(particles.size());
    scalar sumR = 0.;
    forAll(particles, i)
    {
        const mcParticle& p = *particles[i];
        r[i] = mag(p.position() - (p.position()&axis_)*axis_);
        sumR += r[i];
    }
    scalar alpha = particles.size()/sumR;
    forAll(particles, i)
    {
        particles[i]->m() *= alpha*r[i];
    }
}


void Foam::mcParticleCloud::computeCourantNo(mcParticle& p) const
{
    p.Co() = 0.;
    // Due to particle path integration
    const polyMesh& mesh = pMesh();
    const cell& c = mesh.cells()[p.cell()];
    const vector& U = p.Utracking();
    forAll(c, cellFaceI)
    {
        label faceI = c[cellFaceI];
        vector coeff;
        if (mesh.isInternalFace(faceI))
        {
            coeff = CourantCoeffs_[faceI];
        }
        else
        {
            label patchI = mesh.boundaryMesh().whichPatch(faceI);
            if (!CourantCoeffs_.boundaryField()[patchI].size())
            {
                // skip empty boundaries
                continue;
            }
            label i = mesh.boundaryMesh()[patchI].whichFace(faceI);
            coeff = CourantCoeffs_.boundaryField()[patchI][i];
        }
        p.Co() = max(p.Co(), mag(coeff&U));
    }
}


// ************************************************************************* //
