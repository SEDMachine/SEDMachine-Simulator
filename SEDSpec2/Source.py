#!/usr/bin/env python
# 
#  Source.py
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2012-02-02.
#  Copyright 2012 Alexander Rudy. All rights reserved.
#  Version 0.0
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
    from AstroObject.AnalyticSpectra import BlackBodySpectrum,FlatSpectrum,ResampledSpectrum,FLambdaSpectrum
    from AstroObject.AstroSpectra import SpectraObject
    import AstroObject.Utilities as AOU
    import AstroObject.AstroSimulator
except ImportError:
    print "Ensure you have AstroObject installed: please run get-AstroObject.sh"
    raise

import scipy.signal

class SourceCreator(AstroObject.AstroSimulator.Simulator):
    """docstring for SourceCreator"""
    def __init__(self,config=None):
        super(SourceCreator, self).__init__(name="Sources")
        if config == None:
            config = {}
        self.config = config
        self._configureDefaults()
        self.initStages()
    
    def _configureDefaults(self):
        """Sets up the default configuration for the source routines"""
        self.dconfig = {}
        self.dconfig["ExpTime"] = 1200
        
        self.dconfig["plot_format"] = ".pdf"
        
        self.dconfig["Sky"] = {}
        self.dconfig["Sky"]["Massey"] = "MasseySky.fits"
        
        self.dconfig["Configs"] = {}
        self.dconfig["Configs"]["Source"] = "SED.source2.config.yaml"
        
        self.dconfig["Dirs"] = {}
        self.dconfig["Dirs"]["Logs"] = "Logs/"
        self.dconfig["Dirs"]["Partials"] = "Partials/"
        self.dconfig["Dirs"]["Caches"] = "Caches/"
        self.dconfig["Dirs"]["Images"] = "Images/"
        
        
        # Logging Configuration
        self.dconfig["logging"] = {}
        self.dconfig["logging"]["console"] = {}
        self.dconfig["logging"]["console"]["enable"] = True
        self.dconfig["logging"]["console"]["format"] = "......%(message)s"
        self.dconfig["logging"]["console"]["level"] = False
        self.dconfig["logging"]["file"] = {}
        self.dconfig["logging"]["file"]["enable"] = True
        self.dconfig["logging"]["file"]["filename"] = "SEDSource"
        self.dconfig["logging"]["file"]["format"] = "%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s"
        
        AOU.update(self.config,self.dconfig)
        self.config["Configs"]["This"] = self.config["Configs"]["Source"]
        self.log.debug("Set default configuration values.")
    
     
    def initStages(self):
        """Initialize all of the simulator stages"""
        self.registerStage(self.loadSkies,"skydata",description="Loading Sky Data",help="Loading Sky Data")
        self.registerStage(self.setupMoonAttenuation,"moondata",description="Loading Moon Phase Data")
        
    def loadSkies(self):
        """Load in the sky data from files"""
        SKYData = SpectraObject()
        # SKYData.dataClasses = [FLambdaSpectrum]
        SKYData.read("SEDSpec2/MasseySky.fits")
        SKYData.read("SEDSpec2/HansuchikUVES.fits")
        SKYData.read("SEDSpec2/Quimby.fits")
        SKYData.read("SEDSpec2/TurnroseSKY.fits")
        
    def setupMoonAttenuation(self):
        """docstring for setupMoonAttenuation"""
        # Moon phase adjustments. These moon phase attenuation values are for different filter bands.
        # The intermediate wavelengths are accounted for using a polyfit
        # the result is a function which takes phase and wavelength, and outputs an attenuation...

        # See derivation on pg 83 of SED NB 1 (20 July 2011)
        moon_phase = np.array([0., 0.08, 0.16, 0.24, 0.32, 0.40, 0.50])
        moon_g = np.array([2e-17, 2.1e-17, 2.15e-17, 2.3e-17, 5.3e-17, 1.7e-16, 3.2e-16])
        moon_r = np.array([2.3e-17,2.3e-17,2.3e-17,3.3e-17,3.5e-17,8.3e-17,1.3e-16])
        moon_i = np.array([2.8e-17,3.0e-17,3.0e-17,3.3e-17,3.8e-17,7.0e-17,9.0e-17])

        sky_ls = (4868., 6290., 7706., 10000)

        self.moon_funs = []
        for i in xrange(len(moon_phase)):
            gm = moon_g[i]-moon_g[0]
            rm = moon_r[i]-moon_r[0]
            im = moon_i[i]-moon_i[0]
            zm = im

            ff= np.poly1d(np.polyfit(sky_ls, np.array([gm, rm, im, zm]), 2))

            self.moon_funs.append(ff)
        
    
        
        
        
