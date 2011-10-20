// check convergence

Info<< "maxFVResidual = " << maxFVResidual << endl;
Info<< "maxPDFResidual = " << maxPDFResidual << endl;

if (maxFVResidual < FVConvergenceCriterion
    && maxPDFResidual < PDFConvergenceCriterion)
{
    Info<< "reached convergence criterion FV: " << FVConvergenceCriterion
        << "    PDF: " << PDFConvergenceCriterion << endl;
    runTime.writeAndEnd();
    Info<< "latestTime = " << runTime.timeName() << endl;
}

