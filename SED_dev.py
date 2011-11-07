#!/usr/bin/env python
# 
#  SED_dev.py
#  Simulation Software
#  
#  Created by Alexander Rudy on 2011-11-07.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 
# This file is for development with SED Model, and as such, is updated way too frequently.

import numpy as np
import pyfits as pf
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

import SED
from AstroObject.AnalyticSpectra import BlackBodySpectrum, AnalyticSpectrum, FlatSpectrum
import scipy.signal

np.set_printoptions(precision=3,linewidth=100)
lenslet = 2129

System = SED.Model()
Spectrum = BlackBodySpectrum(5000)

System.setup()

System.cache_sed_subimage(lenslet,Spectrum)

System.place_cached_sed(lenslet,"Included Spectrum %d" % lenslet,"Blank")

# System.keep("Included Spectrum %d" % lenslet)
System.write("Experiment.fits",clobber=True)