
    Info<< "\nReading field U\n" << endl;


    // This is a place holder for compressible flows;
    // uniformly 1.0 for now.
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh, 
        dimensionedScalar("rho", dimDensity, 1.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#include "createPhi.H"

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );



    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    volVectorField gradP = fvc::grad(p);

    word mcParticleCloudName("mcParticleCloud");
    args.optionReadIfPresent("cloudName", mcParticleCloudName);

    Info<< "Constructing mcParticleCloud " << mcParticleCloudName << endl;
    IOdictionary particleProperties
    (
       IOobject
        (
            "particleProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    mcParticleCloud mcCloud
    (
     mesh, 
     particleProperties,
     &U,
     &rho,
     &k,
     mcParticleCloudName
    );
