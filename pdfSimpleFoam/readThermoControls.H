bool FVCycle = false;
{
    dictionary thermo= mesh.solutionDict().subDict("thermo");
    label nFVSubCycles = thermo.lookupOrDefault<label>("nFVSubCycles", 1);
    label nPDFSubCycles = thermo.lookupOrDefault<label>("nPDFSubCycles", 1);
    label iCycle = (runTime.timeIndex() - runTime.startTimeIndex() - 1)
        % (nPDFSubCycles + nFVSubCycles);
    FVCycle = iCycle < nFVSubCycles;
}
