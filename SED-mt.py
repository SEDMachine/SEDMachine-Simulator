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

import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool, Value

import SED
from AstroObject.AnalyticSpectra import BlackBodySpectrum
import arpytools.progressbar
import scipy.signal

import math, copy, sys, time, logging, os, argparse

np.set_printoptions(precision=3,linewidth=100)


System = SED.Model()
Spectrum = BlackBodySpectrum(5000)
System.setup()

THREAD_ME = True

lenslets = System.lenslets

DateTime = time.strftime("%Y-%m-%d")
Filename = "Generate-mt-%s" % DateTime



# Setup for Multi-Threading
bar = arpytools.progressbar.ProgressBar()
total = float(len(lenslets))
prog = Value('f',0.0)

def cache_sed(i):
    success = False
    try:
        small, corner = System.get_sub_image_fast(i,Spectrum)
    except SED.SEDLimits:
        SED.LOG.info("Skipped Lenslet %d, Limits Out of Bounds" % i)
    else:
        success = True
        SED.LOG.info("Cached Lenslet %d" % i)
    finally:
        prog.value += 1.0
        bar.render(int((prog.value/total) * 100),"L:%4d" % i)
    if success:
        return (small,corner,i)
    else:
        return False

def cache_sed_cb(results):
    """Takes the result and saves it to the system object"""
    for result in results:
        if hasattr(result,'get'):
            result = result.get(timeout=1)
        if result == False:
            return
        small, corner, lenslet = result
        label = "SUBIMG%d" % lenslet
        System.save(small,label)
        System.frame().metadata=dict(Lenslet=lenslet,Corner=corner)
        SED.LOG.info("Saved Lenslet Cache %d" % lenslet)
    
print("Placing Spectra in %d Lenslets" % len(lenslets))
bar.render(0,"L:%4d" % 0)

# Do the actuall multi-threading

if THREAD_ME:
    pool = Pool()
    pool.map_async(cache_sed,lenslets,1,cache_sed_cb)
    pool.close()
    pool.join()
else:
    results = [cache_sed(i) for i in lenslets]
    cache_sed_cb(results)


print("Rendering Spectra to Full Image...")
bar.render(0,"L:%4d" % 0)
prog = Value('f',0.0)
for i in lenslets:
    try:
        System.place_cached_sed(i,"Included Spectrum %d" % i,"Blank")
        System.show("Blank")
        plt.savefig("Images/Fullimage-%04d.png" % i)
        plt.clf()
        System.show("SUBIMG%d" % i)
        plt.savefig("Images/Subimage-%04d.png" % i)
        plt.clf()
    except SED.SEDLimits:
        SED.LOG.info("Encoutered Unplaced Spectrum %d" % i)
    else:
        SED.LOG.info("Placed Spectrum %d" % i)
    finally:
        prog.value += 1.0
        bar.render(int((prog.value/total) * 100),"L:%4d" % i)
        

System.keep("Blank")
System.select("Blank")
System.crop(System.center[0],System.center[1],1024,1024)
System.write(Filename+".fits",clobber=True)