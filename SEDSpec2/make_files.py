#!/usr/bin/env python
# 
#  make_files.py
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2012-02-02.
#  Copyright 2012 Alexander Rudy. All rights reserved.
#  Version 0.2.0a1
# 

import math
import copy
import sys
import time
import logging
import os
import argparse
import yaml
import collections

import numpy as np
import pyfits as pf
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt

import logging.handlers

import arpytools.progressbar

try:
    from AstroObject.AstroSpectra import SpectraObject
    import AstroObject.Utilities as AOU
    import AstroObject.AstroSimulator
except ImportError:
    print "Ensure you have AstroObject installed: please run get-AstroObject.sh"
    raise

import scipy.signal

# Background Functions
def abmag_to_flambda(AB , lam):
    # Ab magnitude
    # Wavelength in angstrom
    c = 2.9979e18 # angstrom/s
    # return erg/s/cm^2/ang
    return 10**(-(AB + 48.6)/2.5)*c/lam**2

hc = 1.98644521e-8 # erg angstrom

# Moon phase adjustments. These moon phase attenuation values are for different filter bands.
# The intermediate wavelengths are accounted for using a polyfit
# the result is a function which takes phase and wavelength, and outputs an attenuation...

# See derivation on pg 83 of SED NB 1 (20 July 2011)
moon_phase = np.array([0., 0.08, 0.16, 0.24, 0.32, 0.40, 0.50])
moon_g = np.array([2e-17, 2.1e-17, 2.15e-17, 2.3e-17, 5.3e-17, 1.7e-16, 3.2e-16])
moon_r = np.array([2.3e-17,2.3e-17,2.3e-17,3.3e-17,3.5e-17,8.3e-17,1.3e-16])
moon_i = np.array([2.8e-17,3.0e-17,3.0e-17,3.3e-17,3.8e-17,7.0e-17,9.0e-17])

sky_ls = (4868., 6290., 7706., 10000)
moon_funs = []
for i in xrange(len(moon_phase)):
    gm = moon_g[i]-moon_g[0]
    rm = moon_r[i]-moon_r[0]
    im = moon_i[i]-moon_i[0]
    zm = im

    ff= np.poly1d(np.polyfit(sky_ls, np.array([gm, rm, im, zm]), 2))
    
    plt.plot(sky_ls, [gm, rm, im, zm], 'o')
    ll = np.arange(3700, 10000)
    plt.plot(ll, ff(ll))

    moon_funs.append(ff)


# Massey.py -> MasseySky.FITS
from Data.massey import skyab
# Convert some sky spec from AB to fLambda
# This sky spectrum is imported from the massey module, and uses the mag conversion functions above.
skyflam = np.array([skyab[:,0], abmag_to_flambda(skyab[:,1], skyab[:,0])])
MasseySky = SpectraObject(filename="../Data/MasseySky.fits")
MasseySky.save(skyab.T,"SkyAB")
MasseySky.save(skyflam,"SkyFL")
MasseySky.write(clobber=True)

# Hansuchik.py -> HansuchikUVES.FITS
from Data.Hansuchik import uves_sky
HansuchikUVES = SpectraObject(filename="../Data/UVESSky.fits")
HansuchikUVES.save(uves_sky.T,"UVESSky")
HansuchikUVES.write(clobber=True)

# Turnrose.py -> Turnrose.FITS
from Data.Turnrose import skyspec
TurnroseSKY = SpectraObject(filename="../Data/TurnroseSKY.fits")
TurnroseSKY.save(skyspec.T,"TurnroseSKY")
TurnroseSKY.write(clobber=True)

# Quimby.py -> QuimbySky.FITS
from Data.Quimby import quimby_sky
QuimbySKY = SpectraObject(filename="../Data/QuimbySky.fits")
QuimbySKY.save(quimby_sky.T,"QuimbySky")
QuimbySKY.write(clobber=True)

# atmosphere.py -> atmosphere.FITS
from Data.atmosphere import palextinct
atmosphereEXT = SpectraObject(filename="../Data/atmosphere.fits")
atmosphereEXT.save(palextinct.T,"Atmosph")
atmosphereEXT.write(clobber=True)

# palsky_100318.dat -> Palsky.FITS
palskydata = np.genfromtxt("../Data/palsky_100318.dat").T
palSKY = SpectraObject(filename="../Data/PalSky.fits")
palSKY.save(palskydata,"PalSky")
palSKY.write(clobber=True)


