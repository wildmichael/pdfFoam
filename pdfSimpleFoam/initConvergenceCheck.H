// initialize values for convergence checks

    scalar eqnResidual = 1, maxFVResidual = 0, maxPDFResidual = 0;
    scalar FVConvergenceCriterion = 0, PDFConvergenceCriterion = 0;

    simple.readIfPresent("convergence", FVConvergenceCriterion);
    thermo.readIfPresent("convergence", PDFConvergenceCriterion);

