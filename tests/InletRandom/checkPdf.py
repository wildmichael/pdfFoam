#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf, erfc
import scipy.optimize as so
import scipy.stats.distributions as sd
import sys

# config
########
UInSampleMC = np.zeros(1e5) # velocity samples buffer
nInSample = 0             # number of collected samples
npart = 1e6             # number of particles to create
nbins = 100             # number of histogram bins
Umean = 2.0             # mean velocity
urms = 4.0              # velocity variance (rms of fluctuating velocities)
cfl = 1e-1              # particle cfl number

# collect samples using MC
##########################
print "       in"
print "================="
while nInSample < len(UInSampleMC):
  # sample velocities
  velo = (np.random.randn(npart)*urms+Umean)
  # sample positions and make time step
  pos = np.random.rand(npart) + velo*cfl/np.abs(velo).max()
  # find particles with position > 1
  idxIn = np.nonzero(pos > 1.)[0]
  if len(idxIn) and nInSample < len(UInSampleMC):
    # append velocity samples
    if nInSample + len(idxIn) > len(UInSampleMC):
      idxIn = idxIn[:len(UInSampleMC)-nInSample]
    UInSampleMC[nInSample:nInSample+len(idxIn)] = velo[idxIn]
    nInSample += len(idxIn)
  # just print some stats
  print "%7d [%6.2f%%]"%(
      nInSample, 100.*nInSample/len(UInSampleMC))
  sys.stdout.flush()
  # if collected required number of samples, stop
  if nInSample == len(UInSampleMC): break

# collect samples using generator
#################################
import generator
rnd = generator.RandInlet(Umean, 1./(np.sqrt(2.)*urms))
UInSampleGen = np.array([rnd() for i in range(len(UInSampleMC))])

# plotting
##########
Vin = np.linspace(0, Umean+8*urms, 1000)
fig = plt.figure(1)

ax1 = fig.add_subplot(211)
ax1.set_title('inlet (MC)')
ax1.hist(UInSampleMC, bins=nbins, normed=True)
ax1.plot(Vin,rnd._f(Vin,0.))

ax2 = fig.add_subplot(212)
ax2.set_title('inlet (generator)')
ax2.hist(UInSampleGen, bins=nbins,normed=True)
ax2.plot(Vin,rnd._f(Vin,0.))

plt.show()
