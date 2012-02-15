#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  SEDSource.py
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2011-12-03.
#  Copyright 2011 Alexander Rudy. All rights reserved.
#  Version 0.3.0
# 

import math, copy, sys, time, logging, os, argparse, yaml, collections


import math, copy, sys, time, logging, os, argparse

import numpy as np
import pyfits as pf
import scipy as sp

import logging.handlers

import arpytools.progressbar

from multiprocessing import Pool, Value

try:
    from AstroObject.AnalyticSpectra import BlackBodySpectrum,FlatSpectrum,ResampledSpectrum
    from AstroObject.AstroSpectra import SpectraObject
    import AstroObject.Utilities as AOU
    import AstroObject.AstroSimulator
except ImportError:
    print "Ensure you have AstroObject installed: please run get-AstroObject.sh"
    raise

from Utilities import *

import scipy.signal

__version__ = AOU.getVersion(__name__)

class Source(AstroObject.AstroSimulator.Simulator):
    """An SED Machine Source wrapping class"""
    def __init__(self,config):
        super(Source, self).__init__(name="Source")
        self.config = config
        self.debug = self.config["Debug"]
        self.cache = self.config["Cache"]
        self.plot = self.config["Plot"]
        self.defaults = []
        self._configureDefaults()
        self.initStages()
        
    def update(self, d, u):
        """A deep update command for dictionaries.
        This is because the normal dictionary.update() command does not handle nested dictionaries."""
        if u==None:
            return d
        for k, v in u.iteritems():
            if isinstance(v, collections.Mapping):
                r = self.update(d.get(k, {}), v)
                d[k] = r
            else:
                d[k] = u[k]
        return d
    
    def _configureDefaults(self):
        """Sets up the default configuration for the source routines"""
        self.dconfig = {}
        self.dconfig["Type"] = None
        self.dconfig["PreAmp"] = 1.0
        self.dconfig["ExpTime"] = 1200
        
        self.dconfig["plot_format"] = ".pdf"
        
        self.dconfig["Fill"] = {}
        self.dconfig["Fill"]["Type"] = "Flat"
        self.dconfig["Fill"]["PreAmp"] = 1.0
        
        self.dconfig["Object"] = {}
        self.dconfig["Object"]["Type"] = None
        self.dconfig["Object"]["PreAmp"] = 1.0 
        
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
        
        self.defaults += [copy.deepcopy(self.config)]
        self.update(self.config,self.dconfig)
        self.config["Configs"]["This"] = self.config["Configs"]["Source"]
        self.log.debug("Set default configuration values.")
        
    def initStages(self):
        """Initialize simulator stages"""
        self.registerStage(self.setupSource,"sourinit",description="Calculate and resample source")
        self.registerStage(self.setupNoiseandThruput,"thptinit",description="Calculate noise and throughput")
        self.registerMacro("setup","sourinit","thptinit",help="Source initialization functions")
    
    def setup(self):
        """Setup the simulation system"""
        self.commandLine = False
        self.startup()
        self.do("*setup")        
        if self.debug:
            self._plotSpectrum()

        
    
    def _configureFile(self):
        """Setup the configuration, validating the configuration provided"""
        FileName = self.config["Configs"]["Source"]
        try:
            stream = open(FileName,'r')
        except IOError:
            self.log.warning("Configuration File Not Found: %s" % FileName)
        else:
            self.update(self.config,yaml.load(stream))
            stream.close()
            self.log.debug("Loaded Configuration from %s" % FileName)
        finally:
            self.defaults += [copy.deepcopy(self.config)]
    
    def _require(self,element):
        """Require an element"""
        msg = "Configure fails to validate, missing key self.config[\"Source\"][\"%s\"]"
        assert element in self.config, msg % element
    
    def _configureValidate(self):
        """Validation features for the configuration setup"""
        
        if self.config["Type"] == "BlackBody":
            self._require("Temp")
        elif self.config["Type"] == "Flat":
            self._require("Value")
        elif self.config["Type"] == "File":
            self._require("Filename")
        else:
            raise AssertionError("Invalid Type of Source: %s" % self.config["Type"])
        
        
    
    def setupSource(self):
        """A switch function to setup the correct source"""
        self.log.info("Setting up source type '%s'" % self.config["Type"])
        if self.config["Type"] == "BlackBody":
            self._setupBlackbody()
        elif self.config["Type"] == "Flat":
            self._setupFlat()
        elif self.config["Type"] == "File":
            self._setupFile()

    def getSpectrum(self,i):
        """Get spectrum for a particular lenslet"""
        return self.Spectrum
    
    
    def _setupBlackbody(self):
        """Set up a blackbody spectrum"""
        
        self.D_Spectrum = BlackBodySpectrum(self.config["Temp"])
        self.D_Spectrum *= self.config["PreAmp"]
        
    def _setupFlat(self):
        """Sets up a flat spectrum"""
        self.D_Spectrum = FlatSpectrum(self.config["Value"])
        self.D_Spectrum *= self.config["PreAmp"]
        
    
    def _setupFile(self):
        """Sets up a uniform source file based spectrum"""
        WL,Flux = np.genfromtxt(self.config["Filename"]).T
        WL *= 1e-10
        Flux /= np.max(Flux) 
        self.R_Spectrum = ResampledSpectrum(np.array([WL,Flux]),self.config["Filename"])
        self.D_Spectrum = self.R_Spectrum * self.config["PreAmp"]
        
    def setupNoiseandThruput(self):
        """Sets up the model for sky, throughput etc for the isntrument, using Nick's simulation"""
        import SEDSpec.sim
        
        self.log.info("Simulating Noise, Sky and Throughput for System")
        
        lambdas, nsource_photon, shot_noise = SEDSpec.sim.calculate(self.config["ExpTime"],"PI",plot=False,verbose=False)
        
        lambdas = lambdas * 1e-10
        
        l_noise = lambdas[np.isfinite(shot_noise)]
        shot_noise = shot_noise[np.isfinite(shot_noise)]
        l_phot = lambdas[np.isfinite(nsource_photon)]
        nsource_photon = nsource_photon[np.isfinite(nsource_photon)]
        
        
        self.SkySpec = ResampledSpectrum(np.array([l_noise,shot_noise]),"Sky")
        self.SourceFilter = ResampledSpectrum(np.array([l_phot,nsource_photon]),"Source")
        
        self.D_Spectrum *= self.SourceFilter
        self.Spectrum = self.D_Spectrum + self.SkySpec
        
    def _plotSpectrum(self):
        """Plot the spectrum partials"""
        if not self.debug:
            return
        self.log.info("Plotting Source and Intermediate Partials")
        import matplotlib.pyplot as plt
        
        WL = np.linspace(3800,9400,200) * 1e-10
        RE = np.ones(WL.shape) * 100.0
        
        RenderedSpectrum = SpectraObject()
        RenderedSpectrum.save(self.R_Spectrum.data,"Original")
        RenderedSpectrum.save(self.SkySpec.data,"Sky")
        RenderedSpectrum.save(self.SourceFilter.data,"Throughput")
        
        RenderedSpectrum.save(self.D_Spectrum(wavelengths=WL,resolution=RE),"Source")
        
        FL = RenderedSpectrum.data()[1]
        self.log.debug(AOU.npArrayInfo(FL,"RSource"))
        
        
        RenderedSpectrum.save(self.Spectrum(wavelengths=WL,resolution=RE),"Rendered Source")
        
        FL = RenderedSpectrum.data()[1]
        self.log.debug(AOU.npArrayInfo(FL,"RAll"))
        
        filename = "%(directory)sSource-Spectrum%(format)s" % { "directory": self.config["Dirs"]["Partials"], "format": self.config["plot_format"] }
        
        RenderedSpectrum.show()
        plt.title("Source at R=100 constant")
        plt.ylabel("Flux")
        plt.xlabel("Wavelength (m)")
        plt.savefig(filename)
        plt.clf()
        
        filename = "%(directory)sSource-All-Spectra%(format)s" % { "directory": self.config["Dirs"]["Partials"], "format": self.config["plot_format"] }
        RenderedSpectrum.show("Original")
        RenderedSpectrum.show("Sky")
        RenderedSpectrum.show("Throughput")
        RenderedSpectrum.show("Source")
        RenderedSpectrum.show("Rendered Source")
        plt.title("Source at R=100 constant")
        plt.ylabel("Flux")
        plt.xlabel("Wavelength (m)")
        plt.autoscale()
        plt.axis(expandLim(plt.axis()))
        plt.legend()
        plt.savefig(filename)
        plt.clf()
        
        filename = "%(directory)sSource-Spectrum%(format)s" % { "directory": self.config["Dirs"]["Partials"], "format": ".dat" }
        
        np.savetxt(filename,RenderedSpectrum.data("Source"))
        
        filename = "%(directory)sThroughput-Spectrum%(format)s" % { "directory": self.config["Dirs"]["Partials"], "format": ".fits" }
        
        RenderedSpectrum.write(filename,clobber=True)
    
    def dumpConfig(self):
        """Dumps a valid configuration file on top of any old configuration files. This is useful for examining the default configuration fo the system, and providing modifications through this interface."""
        with open(self.config["Configs"]["Source"].rstrip(".yaml")+".dump.yaml",'w') as stream:
            yaml.dump(self.defaults[2]["Source"],stream,default_flow_style=False)


      
