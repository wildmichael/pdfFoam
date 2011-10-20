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
#include "interpolationCellPoint.H"
#include "fixedValueFvPatchField.H"
#include "boundBox.H"
#include "fvc.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "compressible/LES/LESModel/LESModel.H"

// * * * * * * * * * * * * * Local Helper Functions  * * * * * * * * * * * * //

namespace
{
const Foam::compressible::turbulenceModel* getTurbulenceModel
(
    const Foam::objectRegistry& obr
)
{
    if (obr.foundObject<Foam::compressible::turbulenceModel>("RASProperties"))
    {
        return &obr.lookupObject<Foam::compressible::RASModel>("RASProperties");
    }
    else if (obr.foundObject<Foam::compressible::LESModel>("LESProperties"))
    {
        return &obr.lookupObject<Foam::compressible::LESModel>("LESProperties");
    }
    else
    {
        Foam::FatalErrorIn("getTurbulenceModel(const Foam::objectRegistry&)")
            << "No valid model for TKE calculation."
            << Foam::exit(Foam::FatalError);
        return 0;
    }
}


const Foam::dimensionedScalar SMALL_MASS("SMALL_MASS", Foam::dimMass, Foam::SMALL);

const Foam::dimensionedScalar SMALL_VOLUME("SMALL_VOLUME", Foam::dimVolume, Foam::SMALL);

} // anonymous namespace

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineParticleTypeNameAndDebug(mcParticle, 0);
defineTemplateTypeNameAndDebug(Cloud<mcParticle>, 0);
};


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
    volScalarField* rho,
    List<volScalarField*> Phi
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
        U ? *U : mesh.lookupObject<volVectorField>
                     (dict_.lookupOrDefault<word>("U", "U"))
    ),
    rhofv_
    (
        rho ? *rho : const_cast<volScalarField&>(
            mesh.lookupObject<volScalarField>(
                dict_.lookupOrDefault<word>("rho", "rho")))
    ),
    AvgTimeScale_
    (
        dict_.lookupOrAddDefault<scalar>
                ("averageTimeScale", 0.1*runTime_.endTime().value())
    ),
    random_(55555+12345*Pstream::myProcNo()),
    Npc_(dict_.lookupOrAddDefault<label>("particlesPerCell", 30)),
    scalarNames_(0),
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

    ghostCellHash_(256),
    ghostFaceHash_(256),

    pndcPdf_
    (
        IOobject
        (
            "pndCloudPdf",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimDensity,
        mMom_/mesh_.V(),
        rhofv_.boundaryField()
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
        mesh_,
        dimensionedSymmTensor
        (
            "TauCloudPDF", dimVelocity*dimVelocity, symmTensor::zero
        ),
        zeroGradientFvPatchScalarField::typeName
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
        kfv(),
        kfv()().boundaryField()
    ),

    PhicPdf_(Phi),

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
    )
{
    // find scalar fields
    if (PhicPdf_.empty() && dict_.found("scalarFields"))
    {
        dict_.lookup("scalarFields") >> scalarNames_;
        PhicPdf_.setSize(scalarNames_.size());
        forAll(PhicPdf_, PhiI)
        {
            PhicPdf_[PhiI] = &const_cast<volScalarField&>
                (mesh_.lookupObject<volScalarField>(scalarNames_[PhiI]));
        }
    }
    else
    {
        scalarNames_.setSize(PhicPdf_.size());
        forAll(PhicPdf_, PhiI)
        {
            scalarNames_[PhiI] = PhicPdf_[PhiI]->name();
        }
    }
    findGhostLayers();
    checkParticlePropertyDict();
    if (size() > 0) // if particle data found
    {
        mcParticle::readFields(*this);
    }
    else
    {
        Info<< "I am releasing particles initially!" << endl;
        initReleaseParticles();
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
                dimensionedScalar(name, dimless, 0)
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
        updateCloudPDF(0.0);
    }
    else
    {
        FatalErrorIn("Foam::mcParticleCloud::checkMoments()")
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
        PhiMomInstant.set(PhiI, new DimensionedField<scalar, volMesh>
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
            ));
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
    forAllConstIter(mcParticleCloud, *this, pIter)
    {
        const mcParticle& p = pIter();
        vector u = p.UParticle() - p.Updf();
        label cellI = p.cell();

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
    rhofv_.internalField()   = mMom_ / VMomBounded;
    rhofv_.correctBoundaryConditions();

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


// Update particle-cell addressing: for each cell, what are the addresses of the
// particles I host?
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
        else if ( np <= round(Npc_ * cloneAt_) )
        {
            cellPopFlag[celli] = TOOFEW;
        }
        else if (np >= round(Npc_* eliminateAt_) )
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
        if (cellPopFlag[celli] > NORMAL) // either too few or too many particles
        {
            // according to ascending order of mass (if too many particles)
            // or descending order of mass (if too few particle)
            mcParticleList & cepl = cellParticleAddr_[celli];
            if (cepl.size() < 1)
            {
                cepl.setSize(round(PaNIC_[celli]));
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

    //// Debug only: check the list
    //if (debug)
    //{
    //    forAll(cellPopFlag, celli)
    //    {
    //        if (cellPopFlag[celli] == TOOFEW || cellPopFlag[celli] == TOOMANY)
    //        {
    //            Info<< "size is " << cellParticleAddr_[celli].size()
    //                << ", flag is " << cellPopFlag[celli]
    //                << " := " << endl;
    //            mcParticleList & cepl = cellParticleAddr_[celli];

    //            forAll(cepl, pci)
    //            {
    //                if (cepl[pci])
    //                {
    //                    Info<< "ID= " << cepl[pci]->origId() << ", "
    //                        << "m= " << cepl[pci]->m() << endl;
    //                }
    //            }
    //        }
    //    }
    //}
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
    // Pool of partiles to operate on is 2*nx,
    // but liminted by available particles in this cell.
    label nPool = min(2*nx+4, ncur);

    label nKilled = 0;
    scalar massKilled = 0.0;
    for (label particleI=0; particleI < nPool; particleI++)
    {
        mcParticle& p = * cellParticleAddr_[celli][particleI];
        if ( random().scalar01() > 0.5 )
        {
            massKilled += p.m();
            deleteParticle(p);
            cellParticleAddr_[celli][particleI] = 0;
            nKilled++;
        }
        if (nKilled >= nx) break;
    }

    scalar massCompensate = massKilled / (ncur - nKilled);

    // Compensate for the deleted mass
    for (label particleI=0; particleI < ncur; particleI++)
    {
        mcParticle* pPtr =  cellParticleAddr_[celli][particleI];
        if (pPtr)
        {
            pPtr->m() += massCompensate;
        }
    }
    if (debug)
    {
        Info<< "Eliminated " << nKilled
            << " particles in cell " << celli << endl;
    }

}


// Find "ghost cells" (actually the first layer of cells of in/out-flow patch
void Foam::mcParticleCloud::findGhostLayers()
{
    const cellList& cells = mesh_.cells();
    const faceList& faces = mesh_.faces();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    // in/out-flow patches (should read from dictionary)
    const wordList patchNames(dict_.lookup("inoutPatches"));

    label ngPatchs = patchNames.size();
    labelList patchSizes(ngPatchs);

    if (ngPatchs > 0)
    {
        ghostCellLayers_.setSize(ngPatchs);
        ghostFaceLayers_.setSize(ngPatchs);
        ghostCellShifts_.setSize(ngPatchs);
        ghostPatchId_.setSize(ngPatchs);
    }
    else
    {
        return;
    }

    forAll (patchNames, nameI)
    {
        label patchI = patches.findPatchID(patchNames[nameI]);

        if (patchI == -1)
        {
            FatalErrorIn("mcParticleCloud::findGhostLayers()")
                << "Illegal patch " << patchNames[nameI]
                << ". Valid patches are " << patches.name()
                << exit(FatalError);
        }
        label nf = patches[patchI].size();

        patchSizes[nameI] = nf;

        if (nf > 0)
        {
            ghostCellLayers_[nameI].setSize(nf);
            ghostFaceLayers_[nameI].setSize(nf);
            ghostCellShifts_[nameI].setSize(nf);
            ghostPatchId_[nameI] = patchI;
            // Find the ghost cells and faces
            const polyPatch& curPatch = patches[patchI];
            label j = 0;
            forAll(curPatch, facei)
            {
                // find ghost cell
                label faceCelli = curPatch.faceCells()[facei];
                ghostCellLayers_[nameI][j] =  faceCelli;
                //- Add  to hash set
                ghostCellHash_.insert(faceCelli);

                // find opposite face
                label gFaceI = curPatch.start() + facei;
                const cell& ownCell = cells[faceCelli];
                label oppositeFaceI = ownCell.opposingFaceLabel(gFaceI, faces);
                if (oppositeFaceI == -1)
                {
                    FatalErrorIn("mcParticleCloud::findGhostLayers()")
                        << "Face:" << facei << " owner cell:" << faceCelli
                        << " is not a hex?" << abort(FatalError);
                }
                else
                {
                    ghostFaceLayers_[nameI][j] =  oppositeFaceI;
                    ghostCellShifts_[nameI][j] = mesh_.Cf()[oppositeFaceI] - mesh_.Cf().boundaryField()[patchI][facei];
                }

                //- Add face to hash set
                ghostFaceHash_.insert(oppositeFaceI);
                j++;
            }
        }
    }

    // Check the faces found above

    /* forAll(ghostFaceLayers_, patchI)
       forAll(ghostFaceLayers_[patchI], facei)
       {
           Info<<  "face #: " << ghostFaceLayers_[patchI][facei] << endl;
           Info<< faces[ghostFaceLayers_[patchI][facei]].normal(mesh_.points()) << endl;
           Info<< faces[ghostFaceLayers_[patchI][facei]].centre(mesh_.points()) << endl;
       }
       forAll(ghostCellShifts_, patchI)
           forAll(ghostCellShifts_[patchI], celli)
           {Info<< ghostCellShifts_[patchI][celli] << endl;}   */
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mcParticleCloud::evolve()
{

    // lower bound for k
    kcPdf_ = max(kcPdf_, kMin_);

    const volVectorField& gradP = mesh_.lookupObject<const volVectorField>("grad(p)");

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
        fixedValueFvPatchField<scalar>::typeName
    );


    diffRho.internalField() = (pndcPdf_ - rhofv_)/rhofv_;
    diffRho.correctBoundaryConditions();

    volVectorField gradRho = fvc::grad(diffRho) * coeffRhoCorr_;

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

    interpolationCellPoint<scalar> rhoInterp(rhofv_);
    //interpolationCellPointFaceFlux UInterp(U);
    interpolationCellPoint<vector> UInterp(Ufv_);
    interpolationCellPoint<vector> gradPInterp(gradP);
    interpolationCellPoint<scalar> kInterp(kfv());
    interpolationCellPoint<scalar> epsilonInterp(epsilonfv());
    PtrList<interpolationCellPoint<scalar> > PhiInterp(PhicPdf_.size());
    forAll(PhicPdf_, PhiI)
    {
        PhiInterp.set(PhiI,
                      new interpolationCellPoint<scalar>(*PhicPdf_[PhiI]));
    }
    interpolationCellPoint<vector> gradRhoInterp(gradRho);
    interpolationCellPoint<vector> diffUInterp(diffU);
    interpolationCellPoint<scalar> kcPdfInterp(kcPdf_);

    //Impose boundary conditions via particles
    populateGhostCells();

    mcParticle::trackData td(*this, rhoInterp, UInterp, gradPInterp, kInterp,
                             epsilonInterp, PhiInterp, gradRhoInterp,
                             diffUInterp);

    Cloud<mcParticle>::move(td);

    // "Accept" and shift the survived ghost particles
    //  and clear those still in ghost a cell
    purgeGhostParticles();

    scalar existWt = 1.0/(1.0 + (runTime_.deltaT()/AvgTimeScale_).value());
    // Extract statistical averaging to obtain mesh-based quantities
    updateCloudPDF(existWt);
    updateParticlePDF();

    particleNumberControl();

    assertPopulationHealth();
}


// Initialization: populate the FV field with particles
void Foam::mcParticleCloud::initReleaseParticles()
{
    // Populate each cell with Npc_ particles in each cell
    forAll(Ufv_, celli)
    {
        scalar m = mesh_.V()[celli] * rhofv_[celli] / Npc_;
        // TODO fv or pdf?
        vector Updf = Ufv_[celli];
        // TODO shouldn't this be multiplied with 2/3?
        scalar ksqrt = sqrt(kfv()()[celli]);
        vector uscales(ksqrt, ksqrt, ksqrt);
        scalarField Phi(PhicPdf_.size());
        forAll(Phi, PhiI)
        {
            Phi[PhiI] = (*PhicPdf_[PhiI])[celli];
        }

        particleGenInCell(celli, Npc_, m, Updf, uscales, Phi);
    }
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
    boundBox cellbb(pointField(mesh_.points(), mesh_.cellPoints()[celli]), false);

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
                    *this,  position, celli, m, Updf, UParticle, UFap, Phi,
                    runTime_.deltaT().value(), shift, ghost
                );

            addParticle(ptr);
            Npgen ++;
        }

        // until enough particles are generated.
        if(Npgen >= N) break;
    }

    if(Npgen < N)
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


// Enforce in/out flow BCs by populating ghost cells
void Foam::mcParticleCloud::populateGhostCells()
{
    label np = 0;
    label ng = 0;
    forAll(ghostCellLayers_, ghostPatchI)
    {
        forAll(ghostCellLayers_[ghostPatchI], faceCelli)
        {
            ++ ng;
            np += Npc_;

            label  celli = ghostCellLayers_[ghostPatchI][faceCelli];
            scalar m = mesh_.V()[celli] * rhofv_[celli] / Npc_;
            vector Updf = Ufv_[celli];
            // TODO shouldn't this be multiplied with 2/3?
            scalar ksqrt = sqrt(kfv()()[celli]);
            vector uscales(ksqrt, ksqrt, ksqrt);
            label  ghost = 1;
            vector shift = ghostCellShifts_[ghostPatchI][faceCelli];
            // Phi: from patch value (boundary condition)
            scalarField Phi(PhicPdf_.size());
            forAll(Phi, PhiI)
            {
                const volScalarField& f = *PhicPdf_[PhiI];
                label patchId = ghostPatchId_[ghostPatchI];
                Phi[PhiI] = f.boundaryField()[patchId][faceCelli];
            }

            particleGenInCell(celli, Npc_, m, Updf, uscales, Phi, shift, ghost);
        }
    }
    if (debug)
    {
        Info<< np << " particles generated in " << ng << " ghost cells" << endl;
    }

}


void Foam::mcParticleCloud::purgeGhostParticles()
{

    label nDelete = 0;
    label nAdmit = 0;
    forAllIter(mcParticleCloud, *this, pIter)
    {
        mcParticle & p = pIter();
        if (p.ghost() > 0) // This is a ghost
        {
            label celli = p.cell();
            if (ghostCellHash_.found(celli)) // still in ghost cell
            {
                deleteParticle(p);
                ++nDelete;
            }
            else
            {
                // shift and admit this particle as normal member
                p.position() -= p.shift(); // update position
                label newCelli = -1;
                label curCelli = p.cell();
                forAll(mesh_.cellCells()[curCelli], nei)
                {
                    label neiCellId = mesh_.cellCells()[curCelli][nei];
                    if (mesh_.pointInCell(p.position(), neiCellId))
                    {
                        newCelli = neiCellId;
                        break;
                    }
                }
                // Not found in neighbour cells: global search
                if (newCelli < 0)
                {
                    newCelli = mesh_.findCell(p.position());
                }
                if (newCelli < 0)
                {
                    Pout<< "*** Warning in mcParticleCloud::purgeGhostParticles()" << nl
                        << " Shifting of ghost particles caused loss."  << nl
                        << " Possible causes are: " << nl
                        << "  (1) Strange ghost cell shapes;" << nl
                        << "  (2) parallel computing + large time steps" << nl
                        << " Info for the lost particle: " << endl;
                    oneParticleInfo(p);
                    deleteParticle(p);
                    ++nDelete;
                }
                p.cell()  = newCelli;
                p.ghost() = 0;
                p.shift() = vector::zero;
                ++nAdmit;
            }
        }
    }
    Info<< "Ghost particles: " << nDelete << " deleted, "
        << nAdmit << " admitted." << endl;
}

void Foam::mcParticleCloud::info() const
{
    Info<< "Cloud: " << this->name() << nl
        << "    Current number of particles         = "
        << returnReduce(this->size(), sumOp<label>()) << nl;
//      << "    Totally number of particles existed = "
//      << returnReduce(histNp_,      sumOp<label>()) << nl;
}


void Foam::mcParticleCloud::particleInfo() const
{
    forAllConstIter(mcParticleCloud, *this, pIter)
    {
        const mcParticle & p = pIter();
        oneParticleInfo(p);
    }
}


void Foam::mcParticleCloud::oneParticleInfo(const mcParticle& p) const
{
    Pout<< "Particle Id: " << p.origId() << ": "
        << "X     = " << p.position() << ", "
        << "cell  = " << p.cell() << ", "
        << "m     = " << p.m() << nl
        << "U     = " << p.UParticle()  << ", "
        << "Ufv   = " << Ufv_[p.cell()] << ", "
        << "UFap  = " << p.UFap() << ", "
        << "Updf  = " << p.Updf() << nl
        << "Phi   = " << p.Phi() << ", "
        << "ghost = " << p.ghost() << ", "
        << "shift = " << p.shift()
        << endl;
}


void Foam::mcParticleCloud::assertPopulationHealth() const
{
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
            reasons << "\tParticle " << p.origId() << " is a ghost particle" << nl;
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
    if(!cloudHealthy)
    {
        FatalErrorIn("mcParticleCloud::assertPopulationHealth()")
            << "Particle health check not passed:" << nl << reasons.str()
            << exit(FatalError);
    }

}

// ************************************************************************* //
