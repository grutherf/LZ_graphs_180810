#!/usr/bin/env ipython
import libNEST
import numpy as np
import matplotlib.pyplot as plt
import math
#constructor, (particle type,energy,electricfield,xedensity)


NEST  = libNEST.NEST(0,100,1000,2.88,-1) # PT, Energy, E_field, Xe density, dT(us) (use -1 for random uniform dT or -1 in E field for field fringing model)
myDetector = libNEST.Detector()
myDetector.LZSettings()
myDetector.g1 = 0.1 #change g1 to my own value, for fun
NEST.SetDetectorParameters(myDetector)

energies = np.arange(1,50,.01)
ex = []
io = []
photons = []
electrons = []
s1 = []
s2 = []

for energy in energies:

	NEST.SetEnergy(energy)

	NEST.DetectorResponse();

	ex.append(NEST.GetNumExcitons())
	io.append(NEST.GetNumIons())
	photons.append(NEST.GetNumPhotons())
	electrons.append(NEST.GetNumElectrons())
	s1.append(NEST.GetS1())
	s2.append(NEST.GetS2())

plt.xlabel("energy (keV)")
plt.ylabel("S1 (phd)")
plt.scatter(energies,s1)
plt.show()
