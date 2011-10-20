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

SourceFiles
    mcParticleCloudI.H
    mcParticleCloud.C
    mcParticleCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef mcParticleCloud_H
#define mcParticleCloud_H

#include "Cloud.H"
#include "mcParticle.H"
#include "IOdictionary.H"
#include "Random.H"
#include "labelIOField.H"
#include "SortableList.H"
#include "vectorList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

  typedef SortableList<mcParticle*> mcParticleList;
  typedef List <mcParticleList> mcParticleListList;

/*---------------------------------------------------------------------------*\
                           Class mcParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class mcParticleCloud:
    public Cloud<mcParticle>
{

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

    // Private static data

        //- A small mass
        static const dimensionedScalar SMALL_MASS;


    // Private data

        //- The mesh
        const fvMesh& mesh_;
        //- The control dictionary
        dictionary particleProperties_;
        //- The time object
        const Time& runTime_;

        //- The FV velocity field
        const volVectorField& Ufv_;
        //- The FV density field
        const volScalarField& rhofv_;
        //- The FV TKE field
        const volScalarField* kfvPtr_;
        //- The FV mixture fraction field
        volScalarField zfv_;

        //- The averaging time scale
        scalar AvgTimeScale_;

        //- Random number generator for Wiener process (diffusion)
        Random random_;
        //- Desired number of particles per cell (initially)
        label Npc_;
        //- Particle elimination threshold
        scalar eliminateAt_;
        //- Particle cloning threshold
        scalar cloneAt_;
        //- Number of cells
        const label Nc_;
        //- How many particle existed in history (including living ones)
        // only include particles generated in this run.
        scalar histNp_;

        //- List where element @a celli contains a list of particles in cell @celli
        mcParticleListList  cellParticleAddr_;

        // Statistical moments (mass, momentum, energy)

            //- Number of particles in cell
            volScalarField      PaNIC_;
            //- Averaged mass
            volScalarField mMean_;
            //- Averaged volume
            volScalarField VMean_;
            //- Averaged velocity
            volVectorField UMean_;
            //- Averaged mixture fraction
            volScalarField zMean_;
            //- Averaged Reynolds stress
            volSymmTensorField uuMean_;

        //- List of cell-indices for each ghost patch
        labelListList    ghostCellLayers_;
        //- List of labels of the ghost patches
        labelList        ghostPatchId_;
        //- List of face-indices for each ghost patch
        labelListList    ghostFaceLayers_;
        //- List of ghost-cell shifts for each ghost patch
        List<vectorList> ghostCellShifts_;
        //- Hash of ghost cell indices
        HashSet<label>   ghostCellHash_;
        //- Hash of ghost face indices
        HashSet<label>   ghostFaceHash_;

        //- Extracted, time-averaged density
        volScalarField rhocPdf_;
        //- Extracted, time-averaged velocity
        volVectorField UcPdf_;
        //- Extracted, time-averaged mixture fraction
        volScalarField zcPdf_;
        //- Extracted, time-averaged turbulent stress tensor
        volSymmTensorField TaucPdf_;
        //- Extracted, time-averaged TKE
        volScalarField kcPdf_;

        // Model parameters

            //- Density-correction coefficient
            dimensionedScalar coeffRhoCorr_;
            //- Mean velocity relaxation time
            dimensionedScalar URelaxTime_;
            //- Mean TKE relaxation time
            dimensionedScalar kRelaxTime_;
            //- Minimum TKE
            dimensionedScalar kMin_;

        Switch particleNumberControl_;

    // Private Member Functions

        //- Ensure moments are cnosistently read
        void checkMoments();

        //- Update moments and the quantities remembered by particles
        void updateCloudPDF(scalar existWt);

        //- Update Updf & (maybe later, Taupdf) in the particle side.
        void updateParticlePDF();

        // Particle number control actions

            //- Limit cluster/clone thresholds to reasonable values
            void checkParticlePropertyDict();
            //- Perform the particle number control
            void particleNumberControl();
            //- Clone particles in @a celli
            void cloneParticles(label celli);
            //- Eliminate particles in @a celli
            void eliminateParticles(label celli);

        //- Finds the ghost cell layers
        void findGhostLayers();

        //- Enforce in/out flow BCs by populating ghost cells
        void populateGhostCells();

        //- Remove all ghost particles from cloud
        void purgeGhostParticles();

#if 0
        //- Boundary Fields for cloudPDF
        void setCloudPDFBoundary();
        void correctCloudPDFBoundaryConditions();
#endif

        //- Disallow default bitwise copy construct
        mcParticleCloud(const mcParticleCloud&);

        //- Disallow default bitwise assignment
        void operator=(const mcParticleCloud&);


public:

    // Constructors

        //- Construct from components
        // If any of the field pointers happens to be 0, the corresponding
        // field will be looked up from the @a mesh. The names of the fields
        // can be overriden in the dictionary @a dict.
        mcParticleCloud
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const volVectorField* U = 0,
            const volScalarField* rho = 0,
            const volScalarField* k = 0,
            const word& cloudName = "defaultCloud"
        );

    // Member Functions

        //- Access the mesh
        inline const fvMesh& mesh() const;

        //- The TKE FV field
        // NOTE Only safe to call inside evolve().
        inline const volScalarField& kfv() const;
        //- The TKE PDF field
        inline const volScalarField& kcPdf() const;
        //- The TKE relaxation time
        inline const dimensionedScalar& kRelaxTime() const;
        //- The desired number of particles per cell
        inline label Npc();

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
            scalar z,
            scalar rho,
            const vector& shift = vector::zero,
            label  ghost = 0
        );

        //- Evolve the particles
        void evolve();

        //- Report info about name and total number of particles in cloud
        void info() const;
        //- Report info about all particles in cloud
        void particleInfo() const;
        //- Report info about a single particle
        void oneParticleInfo(const mcParticle& p) const;
        //- Asserts that all particles are in domain, none of them is a
        // ghost particle and that its shift is small
        void assertPopulationHealth() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mcParticleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
