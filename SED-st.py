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

from multiprocessing import Pool, Value

import SED
from AstroObject.AnalyticSpectra import BlackBodySpectrum
import arpytools.progressbar
import scipy.signal

np.set_printoptions(precision=3,linewidth=100)


System = SED.Model()
Spectrum = BlackBodySpectrum(5000)
System.setup()

lenslets = np.arange(20) + 2100


# Setup for Multi-Threading
bar = arpytools.progressbar.ProgressBar()
total = float(len(lenslets))
prog = Value('f',0.0)

def cache_sed(i):
    try:
        System.cache_sed_subimage(i,Spectrum)
    except SED.SEDLimits:
        SED.LOG.info("Skipped Lenslet %d, Limits Out of Bounds" % i)
    else:
        SED.LOG.info("Cached Lenslet %d" % i)
    finally:
        prog.value += 1.0
        bar.render(int((prog.value/total) * 100),"L:%d" % i)

print("Placing Spectra in %d Lenslets" % len(lenslets))
bar.render(0,"L:%d" % 0)

for i in lenslets:
    cache_sed(i)

for i in lenslets:
    try:
        System.place_cached_sed(i,"Included Spectrum %d" % i,"Blank")
    except SED.SEDLimits:
        SED.LOG.info("Encoutered Unplaced Spectrum %d" % i)
    else:
        SED.LOG.info("Placed Spectrum %d" % i)

System.keep("Blank")
System.select("Blank")
System.write("Experiment.fits",clobber=True)