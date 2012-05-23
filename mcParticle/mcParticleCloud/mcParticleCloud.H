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

Class
    Foam::mcParticleCloud

Description
    A Cloud of Monte Carlo Particles (mcParticle)

ToDo
    This class uses a second-order time-integration scheme for the particle
    position. However, if particles get reflected at walls, this implementation
    decays to a simple first-order scheme.

    The correct procedure to follow would be to also reflect posititionOld and
    UParticleOld every time a particle hits a wall, such that a particle being
    reflected at a wall is equivalent to that particle entering into the domain
    at the point of reflection. However, positionOld is now outside the domain.
    The procedure to resolve this is as follows:
      1. Find intersection of positionOld + alpha*deltaT*Utracking (alpha being
         a free parameter) with any of the patches on all processors and pick
         the intersection point x that is intersecting from outside with the
         minimum magnitude of alpha.
      2. If alpha < 0 the particle needs to be reflected again at the last wall
         it hit. If after the reflection positionOld is now outside the domain,
         go to 1. Otherwise multiply alpha with -1.
      3. Set positionOld to x, procOld to the processor on which x lies and
         divide Utracking by (1 - alpha) to account for the distance travelled
         outside of the domain during alpha*deltaT.
      3. Proceed as with non-reflected particles.

    This algorithm is very expensive (especially the first step), which is why
    this class currently decays to first-order integration in case of
    reflection.

SourceFiles
    mcParticleCloudI.H
    mcParticleCloud.C
    mcParticleCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef mcParticleCloud_H
#define mcParticleCloud_H

#include "Cloud.H"
#include "mcSolution.H"
#include "DLPtrList.H"
#include "mcParticle.H"
#include "mcBoundary.H"
#include "mcLocalTimeStepping.H"
#include "mcOmegaModel.H"
#include "mcMixingModel.H"
#include "mcPositionCorrection.H"
#include "mcReactionModel.H"
#include "mcVelocityModel.H"
#include "dictionary.H"
#include "Random.H"
#include "labelIOField.H"
#include "SortableList.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "IOLostParticles.H"
#include "richTetPointRef.H"
#include "tetFacePointCellDecomposition.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class mcParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class mcParticleCloud:
    public Cloud<mcParticle>
{

    typedef SortableList<mcParticle*> mcParticleList;
    typedef List <mcParticleList> mcParticleListList;

    //- Particle population status flag
    enum cellPopStatus {EMPTY, NORMAL, TOOFEW, TOOMANY};

    // Private types

        //- scalar comparison functor used for sorting
        class less
        {
            public:

                inline bool operator()
                (
                    mcParticle* one,
                    mcParticle* two
                ) const;
        };

        //- scalar comparison functor used for sorting
        class more
        {
            public:

                inline bool operator()
                (
                    mcParticle* one,
                    mcParticle* two
                ) const;
        };

    // Private data

        //- The mesh
        const fvMesh& mesh_;
        //- The control dictionary
        dictionary thermoDict_;
        //- The solution dictionary
        // Do not access this member directly, use solutionDict() instead.
        mcSolution solutionDict_;
        //- The time object
        const Time& runTime_;
        //- The RAS model object
        const compressible::turbulenceModel& turbModel_;

        //- The FV velocity field
        const volVectorField& Ufv_;

        //- The FV pressure field
        const volScalarField& pfv_;

        //- Random number generator for Wiener process (diffusion)
        Random random_;
        //- List of scalar field names
        wordList scalarNames_;
        //- The scalars to which to apply the mixing model
        labelList mixedScalars_;
        //- Number of cells
        const label Nc_;
        //- How many particle existed in history (including living ones)
        // only include particles generated in this run.
        scalar histNp_;

        //- List where element @a celli contains a list of particles in cell @celli
        mcParticleListList  cellParticleAddr_;

        // Statistical moments (mass, momentum, energy)

            //- Number of particles in cell
            volScalarField PaNIC_;
            //- Averaged 0th statistical moment
            DimensionedField<scalar, volMesh> mMom_;
            //- Averaged 1st statistical moment of inverse-density (volume)
            DimensionedField<scalar, volMesh> VMom_;
            //- Averaged 1st statistical moment of velocity
            DimensionedField<vector, volMesh> UMom_;
            //- Averaged 1st statistical moments of scalar properties
            PtrList<DimensionedField<scalar, volMesh> > PhiMom_;
            //- Averaged 2nd statistical moments of scalar properties
            PtrList<DimensionedField<scalar, volMesh> > PhiRMSMom_;
            //- Averaged 2nd statistical moment of fluctuating velocity
            DimensionedField<symmTensor, volMesh> uuMom_;

        //- List of boundary condition handlers
        PtrList<mcBoundary> boundaryHandlers_;

        //- Owned scalar fields
        //  (i.e. which don't exist already in the registry)
        DLPtrList<volScalarField> ownedScalarFields_;

        //- Extracted, time-averaged density field
        volScalarField& rhocPdf_;
        //- Extracted, instantaneous density field
        volScalarField rhocPdfInst_;
        //- Extracted, instantaneous particle-number density
        volScalarField pndcPdfInst_;
        //- Extracted, time-averaged particle-number density
        volScalarField pndcPdf_;
        //- Extracted, time-averaged velocity
        volVectorField UcPdf_;
        //- Extracted, time-averaged turbulent stress tensor
        volSymmTensorField TaucPdf_;
        //- Extracted, time-averaged TKE
        volScalarField kcPdf_;
        //- Extracted, time-averaged scalar fields
        List<volScalarField*> PhicPdf_;
        //- Extracted, time-averaged RMS scalar fields
        List<volScalarField*> PhiRMScPdf_;

        //- The particle velocity model
        autoPtr<mcVelocityModel> velocityModel_;
        //- The turbulent frequncy model object
        autoPtr<mcOmegaModel> OmegaModel_;
        //- The mixing model object
        autoPtr<mcMixingModel> mixingModel_;
        //- The reaction model object
        autoPtr<mcReactionModel> reactionModel_;
        //- The position correction object
        autoPtr<mcPositionCorrection> positionCorrection_;
        //- The local time stepping object
        autoPtr<mcLocalTimeStepping> localTimeStepping_;

        //- Is this a axi-symmetric case?
        bool isAxiSymmetric_;
        //- Centerline axis for axi-symmetric cases
        vector axis_;
        //- Normal of the centre plane for axi-symmetric cases
        vector centrePlaneNormal_;

        //- Unit-length surface-normal vectors multiplied with the deltaCoeffs
        surfaceVectorField CourantCoeffs_;

        //- Keeps track of the lost particles
        IOLostParticles<mcParticle> lostParticles_;

        //- Scaling factors for the numerical diffusion (~ @c cbrt(mesh.V()))
        scalarList hNum_;

        //- Original, total mass
        scalar m0_;

    // Private Member Functions

        //- Ensure moments are cnosistently read
        void checkMoments();

        //- Update moments and the quantities remembered by particles
        void updateCloudPDF(scalar existWt);

        // Particle number control actions

            //- Perform the particle number control
            void particleNumberControl();
            //- Clone particles in @a celli
            void cloneParticles(label celli);
            //- Eliminate particles in @a celli
            void eliminateParticles(label celli);

        //- Initialize statistical moments
        void initMoments();

        //- Initialize scalar fields
        void initScalarFields();

        //- Initialise boundary condition handlers
        void initBCHandlers();

        //- Compute scaling factors for numerical diffusion
        void initHNum();

        //- Compute and print the mean and maximum particle Courant number
        void printParticleCo();

        //- Generate a @a n random points in cell @a celli
        vectorList randomPointsInCell(label n, label celli);

        //- Disallow default bitwise copy construct
        mcParticleCloud(const mcParticleCloud&);

        //- Disallow default bitwise assignment
        void operator=(const mcParticleCloud&);


public:

    // Constructors

        //- Construct from components
        // If any of the field/object pointers happens to be 0, the
        // corresponding field will be looked up from the @a mesh. The names of
        // the fields can be overriden in the dictionary @a dict.
        mcParticleCloud
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cloudName = "defaultCloud",
            const compressible::turbulenceModel* turbModel = 0,
            const volVectorField* U = 0,
            const volScalarField* p = 0,
            volScalarField* rho = 0
        );

    // Member Functions

        //- Access the mesh
        inline const fvMesh& mesh() const;

        //- Return the thermophysicalProperties dictionary
        inline const dictionary& thermoDict() const;

        //- Return the selected sub-dictionary if the "select" keyword is
        //  given, otherwise return the complete mcSolutionDict
        inline const mcSolution& solutionDict() const;

        //- Read the mcSolution dictionary
        inline virtual bool read();

        //- Read the mcSolution dictionary if modified
        inline virtual bool readIfModified();

        //- Access the turbulence model
        inline const compressible::turbulenceModel& turbulenceModel() const;

        //- The FV velocity field
        inline const volVectorField& Ufv() const;
        //- The TKE FV field
        inline const tmp<volScalarField> kfv() const;
        //- The turbulent dissipation field
        inline const tmp<volScalarField> epsilonfv() const;
        //- The instantaneous PDF density field
        inline const volScalarField& rhocPdfInst() const;
        //- The PDF density field
        inline const volScalarField& rhocPdf() const;
        //- The TKE PDF field
        inline const volScalarField& kcPdf() const;
        //- The scalar PDF fields
        inline const List<volScalarField*>& PhicPdf() const;
        //- The turbulent stress tensor PDF field
        inline const volSymmTensorField& TaucPdf() const;
        //- The particle number (mass) density
        inline const volScalarField& pndcPdf() const;
        //- The instantaneous particle number (mass) density
        inline const volScalarField& pndcPdfInst() const;
        //- The desired number of particles per cell
        inline label Npc();
        //- The names of the scalar fiels
        inline const wordList& scalarNames() const;
        //- List of indices of the scalar fields to which mixing is applied
        inline const labelList& mixedScalars() const;
        //- Whether the case is axi-symmetric
        inline bool isAxiSymmetric() const;
        //- Centerline axis for axi-symmetric cases
        inline const vector& axis() const;
        //- Normal of the centre plane for axi-symmetric cases
        inline const vector& centrePlaneNormal() const;
        //- Coefficients for the particle Courant number
        // (Unit-length surface-normal vectors multiplied with the
        // fvMesh::deltaCoeffs())
        inline const surfaceVectorField& CourantCoeffs() const;

        //- The random number generator
        Random &random() {return random_;}

        //- initial release of particles
        void initReleaseParticles();

        //- Randomly generate N particles in celli, with provided cell-based
        // values and the scale of velocity fluctuation
        void particleGenInCell
        (
            label celli,
            label N,
            scalar m,
            const vector& Updf,
            const vector& uscales,
            const scalarField& Phi,
            const vector& shift = vector::zero,
            label  ghost = 0
        );

        //- Evolve the particles
        // @returns The maximum residual
        scalar evolve();

        //- Handle particles hitting a patch
        template<class TrackData>
        inline void hitPatch
        (
            mcParticle& p,
            TrackData& td,
            label patchI
        );

        //- Report info about name and total number of particles in cloud
        void info() const;
        //- Asserts that all particles are in domain, none of them is a
        // ghost particle and that its shift is small
        void assertPopulationHealth() const;
        //- Post a notification that a particle has been lost
        void notifyLostParticle(const mcParticle& p);

        //- Adjust the mass of generated particles for axi-symmetric cases
        void adjustAxiSymmetricMass(UList<mcParticle*>& particles);

        //- The local time-stepping model
        inline const mcLocalTimeStepping& localTimeStepping() const;

        //- The local time-stepping model
        inline mcLocalTimeStepping& localTimeStepping();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mcParticleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
