    Info<< "Reading field U\n" << endl;
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

    Info<< "Reading field k\n" << endl;
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

    Info<< "Reading field epsilon\n" << endl;
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

    Info<< "Reading field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Creating field pndCloudPDF\n" << endl;
    volScalarField pnd
    (
        IOobject
        (
            "pndCloudPDF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pndCloudPDF", dimDensity, 0.),
        fixedValueFvPatchField<scalar>::typeName
    );

    Info<< "Reading thermophysical properties\n" << endl;
    IOdictionary thermophysicalProperties
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Creating thermo and turbulence models\n" << endl;
    autoPtr<basicRhoThermo> pThermo
    (
        basicRhoThermo::New(mesh)
    );
    basicRhoThermo& thermo = pThermo();

    autoPtr<compressible::RASModel> pTurbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );
    compressible::RASModel& turbulence = pTurbulence();

    Info<< "Creating particle cloud\n" << endl;
    mcParticleCloud cloud
    (
        mesh,
        thermophysicalProperties.subDict("cloudProperties"),
        thermophysicalProperties.lookupOrDefault<word>("cloudName", "mcThermoCloud"),
        &turbulence,
        &U,
        &rho
    );
