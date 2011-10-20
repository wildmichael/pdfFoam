volScalarField rUA = 1.0/UEqn().A();
U = rUA*UEqn().H();
UEqn.clear();

phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rho*rUA, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);

    // Retain the residual from the first iteration
    if (nonOrth == 0)
    {
        eqnResidual = pEqn.solve().initialResidual();
        maxFVResidual = max(eqnResidual, maxFVResidual);
    }
    else
    {
        pEqn.solve();
    }

    if (nonOrth == nNonOrthCorr)
    {
        phi -= pEqn.flux();
    }
}


#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

volVectorField gradP = fvc::grad(p);
U -= rUA*gradP;
U.correctBoundaryConditions();
