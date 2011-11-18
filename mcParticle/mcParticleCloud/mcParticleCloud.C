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
#include "interpolationCellPointFace.H"
#include "fixedValueFvPatchField.H"
#include "boundBox.H"
#include "fvc.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "compressible/LES/LESModel/LESModel.H"
#include "mcProcessorBoundary.H"
#include "gradInterpolationConstantTet.H"

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


//- @todo This is a hack to fix the implementation in OpenFOAM < 2.0.x
template<class T>
void uniqueOrder_FIX
(
    const Foam::UList<T>& lst,
    Foam::labelList& order
)
{
    using namespace Foam;
    sortedOrder(lst, order);

    if (order.size() > 1)
    {
        label n = 0;
        for (label i = 0; i < order.size() - 1; ++i)
        {
            if (lst[order[i]] != lst[order[i+1]])
            {
                order[n++] = order[i];
            }
        }
        // DONT FORGET THE LAST ELEMENT!
        order[n++] = order[order.size()-1];
        order.setSize(n);
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
    return one->m() < two->m();
}


inline bool Foam::mcParticleCloud::more::operator()
(
    mcParticle* one,
    mcParticle* two
) const
{
    return one->m() > two->m();
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
    dict_(dict),
    runTime_(mesh.time()),
    turbModel_
    (
        turbModel ? *turbModel : *getTurbulenceModel(mesh)
    ),
    Ufv_
    (
        U ? *U : mesh_.lookupObject<volVectorField>
                     (dict_.lookupOrDefault<word>("UName", "U"))
    ),
    pfv_
    (
        p ? *p : mesh_.lookupObject<volScalarField>
                     (dict_.lookupOrDefault<word>("pName", "p"))
    ),
    AvgTimeScale_
    (
        dict_.lookupOrAddDefault<scalar>
                ("averageTimeScale", 0.1*runTime_.endTime().value())
    ),
    random_(55555+12345*Pstream::myProcNo()),
    Npc_(dict_.lookupOrAddDefault<label>("particlesPerCell", 30)),
    scalarNames_(0),
    C0_(dict_.lookupOrAddDefault<scalar>("C0", 2.1)),
    C1_(dict_.lookupOrAddDefault<scalar>("C1", 1.0)),
    eliminateAt_(dict_.lookupOrAddDefault<scalar>("eliminateAt", 1.5)),
    cloneAt_(dict_.lookupOrAddDefault<scalar>("cloneAt", 0.67)),
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

    uuMom_
    (
        IOobject
        (
            "uuMoment",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("uuMoment", dimEnergy, symmTensor::zero)
    ),

    ownedScalarFields_(),

    rhocPdf_
    (
        rho ? *rho : const_cast<volScalarField&>(
            mesh_.lookupObject<volScalarField>(
                dict_.lookupOrDefault<word>("rhoName", "rho")))
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

    coeffRhoCorr_
    (
        "cRhoCorr",
        dimTime,
        dict_.lookupOrAddDefault<scalar>("coeffRhoCorrection", 1.0e-7)
    ),
    URelaxTime_
    (
        "URelaxTime",
        dimTime,
        dict_.lookupOrAddDefault<scalar>
                ("URelaxTime", runTime_.deltaT().value()*10.0)
    ),
    kRelaxTime_
    (
        "kRelaxTime",
        dimTime,
        dict_.lookupOrAddDefault<scalar>
                ("kRelaxTime", runTime_.deltaT().value()*10.0)
    ),
    kMin_
    (
        "kMin",
        dimVelocity*dimVelocity,
        dict_.lookupOrAddDefault<scalar>("kMin", 10.0*SMALL)
    ),
    particleNumberControl_
    (
        dict_.lookupOrAddDefault<Switch>("particleNumberControl", true)
    ),
    OmegaModel_(),
    mixingModel_(),
    reactionModel_(),
    isAxiSymmetric_(false),
    CourantCoeffs_
    (
        IOobject
        (
            "CourantCoeffs_",
            mesh_.pointsInstance(),
            mesh_
        ),
        mesh_.surfaceInterpolation::deltaCoeffs() * mesh_.Sf() / mesh_.magSf()
    ),
    lostParticles_(*this)
{
    initScalarFields();

    // Now that the fields exist, create the models
    OmegaModel_ = mcOmegaModel::New(mesh_, dict);
    mixingModel_ = mcMixingModel::New(mesh_, dict);
    reactionModel_ = mcReactionModel::New(mesh_, dict);

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

    initBCHandlers();
    checkParticlePropertyDict();

    // Populate cloud
    if (size() > 0) // if particle data found
    {
        mcParticle::readFields(*this);
    }
    else
    {
        Info<< "I am releasing particles initially!" << endl;
        initReleaseParticles();
    }
    m0_ = 0.;
    forAllConstIter(mcParticleCloud, *this, pIter)
    {
        m0_ += pIter().m();
    }
    // Take care of statistical moments (make sure they are consistent)
    checkMoments();

    // Ensure particles takes the updated PDF values
    updateParticlePDF();
}


// If moments are not read correctly, initialize them.
void Foam::mcParticleCloud::checkMoments()
{
    bool readOk =
        PaNIC_.headerOk() &&
        mMom_.headerOk() &&
        VMom_.headerOk() &&
        UMom_.headerOk() &&
        uuMom_.headerOk();
    // Create moment fields
    PhiMom_.setSize(PhicPdf_.size());
    forAll(PhicPdf_, PhiI)
    {
        // Figure out a field name
        word name = PhicPdf_[PhiI]->name() + "Moment";
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
                dimensionedScalar(name, dimMass, 0)
            ));
        if (!PhiMom_[PhiI].headerOk())
        {
            readOk = false;
        }
    }
    if (readOk)
    {
        Info<< "Moments read correctly." << endl;
    }
    else if (size() > 0)
    {
        Info<< "Moments are missing. Forced re-initialization." << endl;
        initMoments();
    }
    else
    {
        FatalErrorIn("mcParticleCloud::checkMoments()")
            << "Not all moment fields available and no particles present."
            << endl;
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
    PtrList<DimensionedField<scalar, volMesh> > PhiMomInstant(PhiMom_.size());
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
    }
    DimensionedField<symmTensor, volMesh> uuMomInstant
    (
        IOobject
        (
            "uuMomInstant",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("uuMomInstant", dimEnergy, symmTensor::zero)
    );

    PaNIC_ *= 0;

    // Loop through particles to accumulate moments (0, 1, 2 order)
    // as well as particle number
    interpolationCellPointFace<vector> UInterp(Ufv_);
    forAllConstIter(mcParticleCloud, *this, pIter)
    {
        const mcParticle& p = pIter();
        label cellI = p.cell();
        vector U = UInterp.interpolate(p.position(), cellI, p.face());
        vector u = p.UParticle() - U;

        ++PaNIC_[cellI];

        const scalar m = p.m();
        mMomInstant[cellI]  += m;
        VMomInstant[cellI]  += m / p.rho();
        UMomInstant[cellI]  += m * p.UParticle();
        forAll(PhicPdf_, PhiI)
        {
            PhiMomInstant[PhiI][cellI]  += m * p.Phi()[PhiI];
        }
        uuMomInstant[cellI] += m * symm(u * u);
    }

    scalar newWt = 1.0 - existWt;
    // Do time-averaging of moments and compute mean fields
    mMom_  = existWt * mMom_  + newWt * mMomInstant;
    DimensionedField<scalar, volMesh> mMomBounded = max(mMom_, SMALL_MASS);
    pndcPdf_.internalField() = mMom_ / mesh_.V();

    VMom_  = existWt * VMom_  + newWt * VMomInstant;
    DimensionedField<scalar, volMesh> VMomBounded = max(VMom_, SMALL_VOLUME);
    rhocPdf_.internalField()   = mMom_ / VMomBounded;
    rhocPdf_.correctBoundaryConditions();

    UMom_  = existWt * UMom_  + newWt * UMomInstant;
    UcPdf_.internalField()   = UMom_ / mMomBounded;
    UcPdf_.correctBoundaryConditions();

    forAll(PhicPdf_, PhiI)
    {
        PhiMom_[PhiI] = existWt * PhiMom_[PhiI] + newWt * PhiMomInstant[PhiI];
        PhicPdf_[PhiI]->internalField() = PhiMom_[PhiI] / mMomBounded;
        PhicPdf_[PhiI]->correctBoundaryConditions();
    }

    uuMom_ = existWt * uuMom_ + newWt * uuMomInstant;
    TaucPdf_.internalField() = uuMom_/mMomBounded;
    TaucPdf_.correctBoundaryConditions();

    kcPdf_.internalField()   = 0.5 * tr(TaucPdf_.internalField());
    kcPdf_.correctBoundaryConditions();
}


void Foam::mcParticleCloud::updateParticlePDF()
{
    forAllIter(mcParticleCloud, *this, pIter)
    {
        pIter().Updf() = UcPdf_[pIter().cell()];
    }
}


void Foam::mcParticleCloud::checkParticlePropertyDict()
{
    // Cap clone/eliminate threshold with reasonable values
    scalar Dt = runTime_.deltaT().value();

    C0_ = max(0.0, C0_);
    dict_.set("C0", C0_);

    C1_ = max(0.0, min(1.0, C1_));
    dict_.set("C1", C1_);

    eliminateAt_ = max(1.1,  min(eliminateAt_, 2.5));
    dict_.set("eliminateAt", eliminateAt_);

    cloneAt_   = max(0.5,  min(cloneAt_, 0.9));
    dict_.set("cloneAt", cloneAt_);

    kRelaxTime_.value()   = max(2.0*Dt,  kRelaxTime_.value());
    dict_.set("kRelaxTime", kRelaxTime_.value());

    kMin_.value()   = max(1e-3,  min(SMALL, kMin_.value()));
    dict_.set("kMin", kMin_.value());

    URelaxTime_.value()   = max(2.0*Dt,  URelaxTime_.value());
    dict_.set("URelaxTime", URelaxTime_.value());

    Info<< "Sanitized cloudProperties dict:" << dict_ << endl;
}


// Update particle-cell addressing: for each cell, what are the addresses of
// the particles I host?
void Foam::mcParticleCloud::particleNumberControl()
{
    if (!particleNumberControl_) return;

    List<cellPopStatus> cellPopFlag(Nc_, NORMAL);

    forAll(PaNIC_, celli)
    {
        label np = round(PaNIC_[celli]);

        // classify the particle # health condition of each cell
        if (np < 1)
        {
            cellPopFlag[celli] = EMPTY;
        }
        else if ( np < round(Npc_ * cloneAt_) )
        {
            cellPopFlag[celli] = TOOFEW;
        }
        else if (np > round(Npc_* eliminateAt_) )
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
            sort(cellParticleAddr_[celli], more());
            cloneParticles(celli);
        }
        else if ( cellPopFlag[celli] == TOOMANY )
        {
            sort(cellParticleAddr_[celli],  less());
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
    label n = Npc_ - round(PaNIC_[celli]); // no. particle to reproduce
    n = min(round(PaNIC_[celli]), n);

    for (label particleI=0; particleI < n; particleI++)
    {
        mcParticle& p = *(cellParticleAddr_[celli][particleI]);
        // Half my mass
        p.m() /= 2.0;
        // create a new particle like myself
        autoPtr<mcParticle> ptrNew = p.clone();

        addParticle(ptrNew.ptr());
    }
    PaNIC_[celli] += n;

    histNp_ += n;
    if (debug)
    {
        Info<< "Cloned " << n << " particles in cell " << celli << endl;
    }
}


// As name suggests
void Foam::mcParticleCloud::eliminateParticles(label celli)
{
    label ncur = round(PaNIC_[celli]);
    label nx =  ncur - Npc_; // no. particle to eliminate

    // sum of inverse masses
    scalar invM = 0.;
    forAll(cellParticleAddr_[celli], particleI)
    {
        mcParticle& p = *cellParticleAddr_[celli][particleI];
        invM += 1./p.m();
    }

    label nKilled = 0;
    // size of the pool we eliminated particles from
    label nPool = 0;
    scalar massKilled = 0.0;
    // lower bound after which elimination is stopped
    const scalar Nmin = round(cloneAt_ * Npc_);
    forAll(cellParticleAddr_[celli], particleI)
    {
        nPool = particleI;
        mcParticle& p = *cellParticleAddr_[celli][particleI];
        // elimination probability
        scalar Pelim = nx / (p.m() * invM);
        if (random().scalar01() < Pelim)
        {
            ++nKilled;
            massKilled += p.m();
            deleteParticle(p);
            cellParticleAddr_[celli][particleI] = NULL;
        }
        // emergency stop
        if ((ncur - nKilled) < Nmin)
        {
            break;
        }
    }
    ++nPool;
    PaNIC_[celli] -= nKilled;

    // Compensate for the deleted mass
    scalar massCompensate = massKilled / (nPool - nKilled);
    for (label particleI=0; particleI < nPool; particleI++)
    {
        mcParticle* pPtr =  cellParticleAddr_[celli][particleI];
        if (pPtr)
        {
            pPtr->m() += massCompensate;
        }
    }
    if (debug)
    {
        Info<< "Eliminated " << nKilled << " of " << ncur
            << " particles in cell " << celli << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::mcParticleCloud::evolve()
{

    // lower bound for k
    kcPdf_ = max(kcPdf_, kMin_);

    const volVectorField gradP =
        fvc::grad(pfv_ - 2./3.*rhocPdf_*kfv());

    volScalarField diffRho
    (
        IOobject
        (
            "diffRho",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("diffRho", dimless, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    );


    diffRho.internalField() = coeffRhoCorr_*(pndcPdf_ - rhocPdf_)/rhocPdf_;
    diffRho.correctBoundaryConditions();

    volVectorField diffU
    (
        IOobject
        (
            "diffU",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimVelocity/dimTime,
        zeroGradientFvPatchScalarField::typeName
    );

    diffU = (Ufv_ - UcPdf_) / URelaxTime_;
    diffU.correctBoundaryConditions();

    interpolationCellPointFace<scalar> rhoInterp(rhocPdf_);
    //interpolationCellPointFaceFlux UInterp(U);
    interpolationCellPointFace<vector> UInterp(Ufv_);
    interpolationCellPointFace<vector> gradPInterp(gradP);
    interpolationCellPointFace<scalar> kInterp(kfv());
    gradInterpolationConstantTet<scalar> gradRhoInterp(diffRho);
    interpolationCellPointFace<vector> diffUInterp(diffU);
    interpolationCellPointFace<scalar> kcPdfInterp(kcPdf_);

    // Correct boundary conditions
    forAll(boundaryHandlers_, boundaryI)
    {
        boundaryHandlers_[boundaryI].correct(*this, false);
    }

    OmegaModel_().correct(*this);
    mixingModel_().correct(*this);
    reactionModel_().correct(*this);

    scalar existWt = 1.0/(1.0 + (runTime_.deltaT()/AvgTimeScale_).value());

    mcParticle::trackData td
    (
        *this,
        rhoInterp,
        UInterp,
        gradPInterp,
        kInterp,
        gradRhoInterp,
        diffUInterp
    );

    forAllIter(mcParticleCloud, *this, pIter)
    {
        pIter().nSteps() = 0;
    }

    Cloud<mcParticle>::move(td);

    // Correct boundary conditions
    forAll(boundaryHandlers_, boundaryI)
    {
        boundaryHandlers_[boundaryI].correct(*this, true);
    }

    // Extract statistical averaging to obtain mesh-based quantities
    updateCloudPDF(existWt);
    updateParticlePDF();
    printParticleCo();

    particleNumberControl();

    if (debug)
    {
        assertPopulationHealth();
    }

    scalar rhoRes = 0./*,
           URes   = 0.,
           kRes   = 0.*/;

    forAll(rhocPdf_, cellI)
    {
        rhoRes += fabs((pndcPdf_[cellI]-rhocPdf_[cellI])/rhocPdf_[cellI]);
#if 0
        URes += mag
        (
            cmptDivide
            (
                UcPdf_[cellI] - Ufv_[cellI],
                stabilise(Ufv_[cellI], 100*SMALL)
            )
        );
        scalar k = kfv()()[cellI];
        kRes += fabs((kcPdf_[cellI] - k)/k);
#endif
    }
    reduce(rhoRes, sumOp<scalar>());
#if 0
    reduce(URes, sumOp<scalar>());
    reduce(kRes, sumOp<scalar>());
#endif
    Info<< "Cloud " << name() << " residuals: rho = " << rhoRes
        /*<< ", U = " << URes << ", k = " << kRes */<< endl;

    label nLostPart = lostParticles_.size();
    reduce(nLostPart, sumOp<label>());
    Info<< "Cloud " << name() << " number of lost particles: "
        << nLostPart << nl;
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
    deltaU.write();
    Info<< "DEBUG: maximum deltaU = " << gMax(deltaU) << endl;
#endif
    // Mass after the evolution done
    scalar m1 = 0;
    forAllConstIter(mcParticleCloud, *this, pIter)
    {
        m1 += pIter().m();
    }
    // Difference of the masses
    scalar diffM = m1-m0_;
    Info<< "    m0 is "<< m0_ << endl;
    Info<< "    m1 is "<< m1 << endl;
    Info<< "    The difference of the mass is "<< diffM << endl;
    Info<< "    The relative difference of the mass (in percent) is "
        << diffM/m0_*100  << endl;
    return rhoRes/*+URes+kRes*/;
}


// Initialization: populate the FV field with particles
void Foam::mcParticleCloud::initReleaseParticles()
{
    // Populate each cell with Npc_ particles in each cell
    forAll(Ufv_, celli)
    {
        scalar m = mesh_.V()[celli] * rhocPdf_[celli] / Npc_;
        // TODO fv or pdf?
        vector Updf = Ufv_[celli];
        scalar urms = sqrt(2./3.*kfv()()[celli]);
        vector uscales(urms, urms, urms);
        scalarField Phi(PhicPdf_.size());
        forAll(Phi, PhiI)
        {
            Phi[PhiI] = (*PhicPdf_[PhiI])[celli];
        }

        particleGenInCell(celli, Npc_, m, Updf, uscales, Phi);
    }
    OmegaModel_().correct(*this);
    reactionModel_().correct(*this);
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

    label Npgen = 0;
    for (int i = 0; i < 100 * N; i++)
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

        // Initially put $N particle per cell
        if (mesh_.pointInCell(position, celli))
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
            vector UFap = Ufv_[celli];

            mcParticle* ptr = new mcParticle
            (
                *this,
                position,
                celli,
                m,
                Updf,
                UParticle,
                UFap,
                Phi,
                shift,
                ghost
            );

            addParticle(ptr);
            ++Npgen;
        }

        // until enough particles are generated.
        if (Npgen >= N) break;
    }

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
    mMom_  = mesh_.V() * rhocPdf_;
    pndcPdf_.internalField() = rhocPdf_.internalField();

    VMom_  = mMom_ / rhocPdf_;

    UMom_  = mMom_ * Ufv_;
    UcPdf_.internalField()   = Ufv_.internalField();
    UcPdf_.correctBoundaryConditions();

    forAll(PhicPdf_, PhiI)
    {
        PhiMom_[PhiI] = mMom_ * (*PhicPdf_[PhiI]);
    }

    uuMom_ = mMom_ * turbulenceModel().R()();

    kcPdf_.internalField()   = turbulenceModel().k()().internalField();
    kcPdf_.correctBoundaryConditions();
}


void Foam::mcParticleCloud::initScalarFields()
{
    // Find scalar fields
    if (dict_.found("scalarFields"))
    {
        dict_.lookup("scalarFields") >> scalarNames_;
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
                << dict_.lookupEntry("scalarFields", false, false).name()
                << " contains duplicate entries.\n"
                << exit(FatalError);
        }
    }
    label nScalarFields = scalarNames_.size();
    PhicPdf_.setSize(nScalarFields);
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
    }

    // Find labels of mixed scalars
    scalarNames_.setSize(nScalarFields);
    wordList mixedScalarNames;
    if (dict_.found("mixedScalars"))
    {
        dict_.lookup("mixedScalars") >> mixedScalarNames;
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
                << dict_.lookupEntry("mixedScalars", false, false).name()
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
}


void Foam::mcParticleCloud::initBCHandlers()
{
    boundaryHandlers_.clear();
    boundaryHandlers_.setSize(mesh_.boundaryMesh().size());
    const dictionary& bd = dict_.subDict("boundaryHandlers");
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        autoPtr<mcBoundary> bcHandler;
        if (isA<processorPolyPatch>(pp))
        {
            bcHandler.reset(new mcProcessorBoundary
            (
                mesh_,
                patchI
            ));
        }
        else
        {
            bcHandler = mcBoundary::New
            (
                mesh_,
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


void Foam::mcParticleCloud::printParticleCo()
{
    scalar meanPartCoNum = 0.0;
    scalar maxPartCoNum = 0.0;

    if (mesh_.nInternalFaces() && size())
    {
        forAllConstIter(mcParticleCloud, *this, pIter)
        {
            meanPartCoNum += pIter().Co();
            maxPartCoNum = max(maxPartCoNum, pIter().Co());
        }
        label nPart = size();
        reduce(nPart, sumOp<label>());
        reduce(meanPartCoNum, sumOp<scalar>());
        meanPartCoNum /= nPart;
        reduce(maxPartCoNum, maxOp<scalar>());
    }

    Info<< "Particle Courant Number in cloud " << name()
        << " mean: " << meanPartCoNum
        << " max: " << maxPartCoNum << endl;
}


void Foam::mcParticleCloud::notifyLostParticle(const Foam::mcParticle& p)
{
    // TODO Distribute mass to other particles in cell
    lostParticles_.add(p);
}

// ************************************************************************* //
