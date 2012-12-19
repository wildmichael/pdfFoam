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

Application
    mcInletRandomTest

Description
    Runs mcInletRandom and dumps the generated numbers to stdout

\*---------------------------------------------------------------------------*/

#include "mcInletRandom.H"
#include "Random.H"
#include "error.H"
#include "labelList.H"
#include "OFstream.H"
#include "cpuTime.H"
#include <cstdio>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

    const scalar Umean = 2.,         // mean velocity (mean(U))
                 urms = 4.,          // RMS of fluctuations (std(U))
                 hMin = 0.0,         // lower cutoff of histogram
                 hMax = 20.0;        // upper cutoff of histogram
    const label  Nbins = 100,        // number of histogram bins
                 Nsamples = 1000000; // number of samples to draw
    const scalar w = (hMax - hMin) / Nbins; // histogram bin width

    Random rnd(0);
    wordList types(2);
    types[0] = "sampling";
    types[1] = "inversion";

    labelListList hist(types.size());

    forAll(types, typeI)
    {
        word t = types[typeI];
        dictionary d;
        d.add("type", t);
        autoPtr<mcInletRandom> irnd = mcInletRandom::New(rnd, Umean, urms, d);

        // Burning in
        for (label i=0; i<100000; ++i)
        {
            scalar v = rnd.scalar01();
            (void)v; // silence warning about unuse variable v
        }

        hist[typeI].setSize(Nbins);       // histogram counts
        for (label i=0; i<Nbins; ++i) hist[typeI][i] = 0;
        // collect samples (only count those in histogram range)
        cpuTime c;
        for (label i=0; i<Nsamples; ++i)
        {
            scalar v = irnd().value();
            if (v < hMin || v > hMax)
            {
                continue;
            }
            // find bin into which this sample belongs
            for (label j=0; j<Nbins; ++j)
            {
                if (v < (w*(j+1)+hMin))
                {
                    ++hist[typeI][j];
                    break;
                }
            }
        }
        Info<< "Generating " << label(Nsamples) << " random numbers with " << types[typeI]
            << " took " << c.elapsedCpuTime() << " seconds.\n";
    }

    // output histogram data (bin center and normalized counts)
    OFstream outf("histogram.dat");
    outf<< "# center";
    forAll(types, typeI)
    {
        outf<< tab << types[typeI];
    }
    outf<< nl;
    for (label i=0; i<Nbins; ++i)
    {
        scalar c = i*w+0.5*w+hMin;
        outf<< c;
        forAll(types, typeI)
        {
           outf<< tab << scalar(hist[typeI][i])/(w*Nsamples);
        }
       outf<< nl;
    }
    return 0;
}

// ************************************************************************* //
