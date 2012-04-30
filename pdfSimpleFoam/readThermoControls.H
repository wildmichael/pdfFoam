bool FVCycle = false;
{
    static label subCyclesCount = 0;
    dictionary cd = runTime.controlDict();
    label nFVSubCycles = readLabel(cd.lookup("nFVSubCycles"));
    label nTotalCycles = nFVSubCycles + readLabel(cd.lookup("nPDFSubCycles"));
    if (subCyclesCount >= nTotalCycles)
    {
        subCyclesCount = 0;
    }
    FVCycle = subCyclesCount < nFVSubCycles;
    ++subCyclesCount;
}
