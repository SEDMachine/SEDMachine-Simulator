#!/usr/bin/env python
# 
#  BlackBodyImage.py
#  Simulation Software
#  
#  Created by Alexander Rudy on 2011-10-21.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 

from SEDImage import *

bbody = BlackBodySpectrum(5000)
bbody1 = BlackBodySpectrum(100)

test = SEDImage()

test.loadOpticsData("xy_13oct2011_v40.TXT","Rtest.dat")
test.generate_blank()

points, wl, deltawl = test.get_wavelengths(2129)

radiance = bbody(wl*1e-6)
flux = radiance[1,:-1] * deltawl

plt.subplot(2,1,1)
plt.plot(wl,radiance[1,:],'b.')
plt.title("Blackbody Spectrum")
plt.ylabel("Spectral Radiance (Radiance per unit wavelength, SI Units W * sr-1 * m-2 * m-1 )")
plt.xlabel("Wavelegnth")
plt.subplot(2,1,2)
plt.plot(wl[:-1],flux,'g.')
plt.ylabel("Radiance (Radiance, SI Units W * sr-1 * m-2 )")
plt.xlabel("Wavelegnth")
plt.savefig("RF.png")
