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
    from AstroObject.AnalyticSpectra import BlackBodySpectrum,FlatSpectrum,ResampledSpectrum,FLambdaSpectrum,InterpolatedSpectrum
    from AstroObject.AstroSpectra import SpectraObject
    import AstroObject.Utilities as AOU
    import AstroObject.AstroSimulator
except ImportError:
    print "Ensure you have AstroObject installed: please run get-AstroObject.sh"
    raise

import scipy.signal

class SourceCreator(AstroObject.AstroSimulator.Simulator):
    """docstring for SourceCreator"""
    def __init__(self):
        super(SourceCreator, self).__init__(name="Sources")
        self.config.merge(self.defaults)
        self.initStages()
    
    defaults = {
        "ExpTime" : 1200,
        "plot_format" : ".pdf",
        "Sky" : {
            "Massey" : "MasseySky.fits",  
        },
        "Moon" : {
            "Phase" : 1
        },
        "Thpt" : "SEDSpec2/Data/thpt.npy",
        "QE" : "prism_pi",
        "Configs" : {
            "Source" : "SED.source2.config.yaml"
        },
        "Dirs": {
            "Logs": "Logs",
            "Partials" : "Partials",
            "Caches" : "Caches",
            "Images" : "Images",
        },
        'logging': {
            'console': {
                'level': logging.INFO, 
                'enable': True, 
                'format': '%(levelname)-8s... %(message)s'
            },
            'growl' : {
                'name' : "SED Machine Simulator",
            },
            'file': {
                'format': '%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s',
                'enable': True, 
                'filename': 'SEDMachine'
            },
        },
    }
     
    def initStages(self):
        """Initialize all of the simulator stages"""
        self.registerStage(self.physical_constants,"constants",description="Setting Physical Constants",help=False)
        self.registerStage(self.load_skies,"sky-data",description="Loading Sky Data",help="Loading Sky Data")
        self.registerStage(self.setup_moon,"moon-data",description="Loading Moon Phase Data",dependencies=["constants"])
        self.registerStage(None,"sky-func",description="DUMMY Setting Sky Function")
        self.registerStage(self.use_turnrose,"turnrose",description="Setting Turnrose Sky Function",help="Use Turnrose Sky spectrum",dependencies=["sky-data","moon-data","constants"],replaces=["sky-func"])
        self.registerStage(self.setup_thpt,"setup-thruput",description="Making Thruput Functions",help=False)
        self.registerStage(self.include_moon,"moon",description="Applying Moon Light",help=False,dependencies=["moon-data","sky-func"])
        self.registerStage(self.make_sky,"make-sky",description="Make final sky spectrum",help=False,dependencies=["setup-thruput"])
        
        
    def physical_constants(self):
        """Setup physical constants"""
        self.const = {}
        self.const["hc"] = 1.98644521e-8 # erg angstrom

        
    def load_skies(self):
        """Load in the sky data from files"""
        self.SKYData = SpectraObject()
        # SKYData.dataClasses = [FLambdaSpectrum]
        self.SKYData.read("SEDSpec2/MasseySky.fits")
        self.SKYData.read("SEDSpec2/HansuchikUVES.fits")
        self.SKYData.read("SEDSpec2/Quimby.fits")
        self.SKYData.read("SEDSpec2/TurnroseSKY.fits")
        
    def setup_moon(self):
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

        self.moon_specs = []
        for i in xrange(len(moon_phase)):
            gm = moon_g[i]-moon_g[0]
            rm = moon_r[i]-moon_r[0]
            im = moon_i[i]-moon_i[0]
            zm = im
            fluxes = np.array([gm, rm, im, zm])/self.const["hc"]/sky_ls
            moon_spec = InterpolatedSpectrum(np.array([sky_ls,fluxes]),"Moon Phase %s" % i)
            
            self.moon_specs.append(moon_spec)
        
    def use_turnrose(self):
        """docstring for setupTurnrose"""
        WL,FL = self.SKYData.data("TurnroseSKY")
        FL *= 1e-18 * 3 # NICK! I NEED THIS EXPLAINED!
        FL /= self.const["hc"] / WL
        self.SkySpectrum = FLambdaSpectrum(np.array([WL*1e-10,FL]),"SkySpectrum")
        
        
    def include_moon(self):
        """docstring for include_moon"""
        self.SkySpectrum += self.moon_specs[self.config["Moon"]["Phase"]]
        
        
        
    def setup_thpt(self):
        """Sets up thruputs"""
        thpts = np.load(self.config["Thpt"])[0]
        self.qe = {}
        self.qe["prism_pi"] = InterpolatedSpectrum(np.array([thpts["lambda"], thpts["thpt-prism-PI"]]),"PI Prism")
        self.qe["prism_andor"] = InterpolatedSpectrum(np.array([thpts["lambda"], thpts["thpt-prism-Andor"]]),"Andor Prism")
        self.qe["grating"] = InterpolatedSpectrum(np.array([thpts["lambda"], thpts["thpt-grating"]]),"Grating")
        
    def make_sky(self):
        """Make a sky spectrum"""
        self.SkySpectrum *= self.qe[self.config["QE"]]
        
    
        
def run():
    SIM = SourceCreator()
    SIM.run()
    
if __name__ == '__main__':
    run()        
        
