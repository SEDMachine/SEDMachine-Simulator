#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  SEDSource.py
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2011-12-03.
#  Copyright 2011 Alexander Rudy. All rights reserved.
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
except ImportError:
    print "Ensure you have AstroObject installed: please run get-AstroObject.sh"
    raise

from Utilities import *

import scipy.signal

__version__ = open(os.path.abspath(os.path.join(os.path.dirname(__file__),"VERSION")),'r').read()

class Source(object):
    """An SED Machine Source wrapping class"""
    def __init__(self,config):
        super(Source, self).__init__()
        self.config = config
        self.debug = self.config["System"]["Debug"]
        self.cache = self.config["System"]["Cache"]
        self.plot = self.config["System"]["Plot"]
        self.defaults = []
        self.initLog()
    
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
        self.config["Source"] = {}
        self.config["Source"]["Type"] = None
        self.config["Source"]["PreAmp"] = 1.0
        self.config["Source"]["ExpTime"] = 1200
        
        self.config["Source"]["plot_format"] = ".pdf"
        
        self.config["Source"]["Fill"] = {}
        self.config["Source"]["Fill"]["Type"] = "Flat"
        self.config["Source"]["Fill"]["PreAmp"] = 1.0
        
        self.config["Source"]["Object"] = {}
        self.config["Source"]["Object"]["Type"] = None
        self.config["Source"]["Object"]["PreAmp"] = 1.0 
        
        # Logging Configuration
        self.config["Source"]["logging"] = {}
        self.config["Source"]["logging"]["console"] = {}
        self.config["Source"]["logging"]["console"]["enable"] = True
        self.config["Source"]["logging"]["console"]["format"] = "......%(message)s"
        self.config["Source"]["logging"]["console"]["level"] = False
        self.config["Source"]["logging"]["file"] = {}
        self.config["Source"]["logging"]["file"]["enable"] = True
        self.config["Source"]["logging"]["file"]["filename"] = "SEDSource"
        self.config["Source"]["logging"]["file"]["format"] = "%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s"
        
        self.defaults += [copy.deepcopy(self.config)]
        self.log.debug("Set default configuration values.")
        
    def initLog(self):
        """Initializes the system logger. This logger starts with only a buffer, no actual logging output. The buffer is used to hold log messages before a logging output location has been specified."""
        self.log = logging.getLogger(__name__)
        self.log.setLevel(logging.DEBUG)
        logging.captureWarnings(True)
        self.filebuffer = logging.handlers.MemoryHandler(1e6) 
        #Our handler will only handle 1-million messages... lets not go that far
        self.consolebuffer = logging.handlers.MemoryHandler(1e6)
        self.consolebuffer.setLevel(logging.INFO)
        self.log.addHandler(self.filebuffer)
        self.log.addHandler(self.consolebuffer)
        
        self.log.info("--------------------------------------")
        self.log.info("Welcome to the SEDMachine Source model")
        self.log.debug("Version %s" % __version__ )
    
    def setup(self):
        """Setup the simulation system"""
        
        self._configureDefaults()
        self._configureSystem()
        self._configureFile()
        
        self.setupLog()
        
        self.setupSource()
        
        self.setupNoiseandThruput()
        
        if self.debug:
            self._plotSpectrum()
    
    def _configureSystem(self):
        """Configure the source based on the system configuration of the soure"""
        
        self.update(self.config["Source"],self.config["System"]["Source"])
        
        self.defaults += [copy.deepcopy(self.config)]
        
    
    def _configureFile(self):
        """Setup the configuration, validating the configuration provided"""
        FileName = self.config["System"]["Configs"]["Source"]
        try:
            stream = open(FileName,'r')
        except IOError:
            self.log.warning("Configuration File Not Found: %s" % FileName)
        else:
            self.update(self.config["Source"],yaml.load(stream))
            stream.close()
            self.log.debug("Loaded Configuration from %s" % FileName)
        finally:
            self.defaults += [copy.deepcopy(self.config)]
    
    def _require(self,element):
        """Require an element"""
        msg = "Configure fails to validate, missing key self.config[\"Source\"][\"%s\"]"
        assert element in self.config["Source"], msg % element
    
    def _configureValidate(self):
        """Validation features for the configuration setup"""
        
        if self.config["Source"]["Type"] == "BlackBody":
            self._require("Temp")
        elif self.config["Source"]["Type"] == "Flat":
            self._require("Value")
        elif self.config["Source"]["Type"] == "File":
            self._require("Filename")
        else:
            raise AssertionError("Invalid Type of Source: %s" % self.config["Source"]["Type"])
        
    
    def setupLog(self):
        """Setup Logging Functions for the SEDMachine Model.

        This configures the logging system, including a possible console and file log. It then reads the logging buffer into the logfile.
        """

        self.console = logging.StreamHandler()
        consoleFormat = self.config["Source"]["logging"]["console"]["format"]
        if self.config["Source"]["logging"]["console"]["level"]:
            self.console.setLevel(self.config["Source"]["logging"]["console"]["level"])
        elif self.debug:
            self.console.setLevel(logging.DEBUG)
        else:
            self.console.setLevel(logging.ERROR)
        consoleFormatter = logging.Formatter(consoleFormat)
        self.console.setFormatter(consoleFormatter)
        if self.config["Source"]["logging"]["console"]["enable"]:
            self.log.addHandler(self.console)
            self.consolebuffer.setTarget(self.console)
        self.consolebuffer.close()
        self.log.removeHandler(self.consolebuffer)
        
        self.logfile = None
        # Only set up the file log handler if we can actually access the folder
        if os.access(self.config["System"]["Dirs"]["Logs"],os.F_OK) and self.config["Source"]["logging"]["file"]["enable"]:
            filename = self.config["System"]["Dirs"]["Logs"] + self.config["Source"]["logging"]["file"]["filename"]+".log"
            self.logfile = logging.handlers.TimedRotatingFileHandler(filename,when='midnight')
            self.logfile.setLevel(logging.DEBUG)
            fileformatter = logging.Formatter(self.config["Source"]["logging"]["file"]["format"],datefmt="%Y-%m-%d-%H:%M:%S")
            self.logfile.setFormatter(fileformatter)
            self.log.addHandler(self.logfile)
            # Finally, we should flush the old buffers
            self.filebuffer.setTarget(self.logfile)
        
        self.filebuffer.close()
        self.log.removeHandler(self.filebuffer)
        self.log.debug("Configured Logging")
        
    
    def setupSource(self):
        """A switch function to setup the correct source"""
        self.log.info("Setting up source type '%s'" % self.config["Source"]["Type"])
        if self.config["Source"]["Type"] == "BlackBody":
            self._setupBlackbody()
        elif self.config["Source"]["Type"] == "Flat":
            self._setupFlat()
        elif self.config["Source"]["Type"] == "File":
            self._setupFile()

    def getSpectrum(self,i):
        """Get spectrum for a particular lenslet"""
        return self.Spectrum
    
    
    def _setupBlackbody(self):
        """Set up a blackbody spectrum"""
        
        self.Spectrum = BlackBodySpectrum(self.config["Source"]["Temp"])
        self.Spectrum *= self.config["Source"]["PreAmp"]
        
    def _setupFlat(self):
        """Sets up a flat spectrum"""
        self.Spectrum = FlatSpectrum(self.config["Source"]["Value"])
        self.Spectrum *= self.config["Source"]["PreAmp"]
        
    
    def _setupFile(self):
        """Sets up a uniform source file based spectrum"""
        WL,Flux = np.genfromtxt(self.config["Source"]["Filename"]).T
        WL *= 1e-10
        Flux /= np.max(Flux) 
        self.R_Spectrum = ResampledSpectrum(np.array([WL,Flux]),self.config["Source"]["Filename"])
        self.D_Spectrum = self.R_Spectrum * self.config["Source"]["PreAmp"]
        
    def setupNoiseandThruput(self):
        """Sets up the model for sky, throughput etc for the isntrument, using Nick's simulation"""
        import SEDSpec.sim
        
        self.log.info("Simulating Noise, Sky and Throughput for System")
        
        lambdas, nsource_photon, shot_noise = SEDSpec.sim.calculate(self.config["Source"]["ExpTime"],"PI",plot=False,verbose=False)
        
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
        
        filename = "%(directory)sSource-Spectrum%(format)s" % { "directory": self.config["System"]["Dirs"]["Partials"], "format": self.config["Source"]["plot_format"] }
        
        RenderedSpectrum.show()
        plt.title("Source at R=100 constant")
        plt.ylabel("Flux")
        plt.xlabel("Wavelength (m)")
        plt.savefig(filename)
        plt.clf()
        
        filename = "%(directory)sSource-All-Spectra%(format)s" % { "directory": self.config["System"]["Dirs"]["Partials"], "format": self.config["Source"]["plot_format"] }
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
        
        filename = "%(directory)sSource-Spectrum%(format)s" % { "directory": self.config["System"]["Dirs"]["Partials"], "format": ".dat" }
        
        np.savetxt(filename,RenderedSpectrum.data("Source"))
        
        filename = "%(directory)sThroughput-Spectrum%(format)s" % { "directory": self.config["System"]["Dirs"]["Partials"], "format": ".fits" }
        
        RenderedSpectrum.write(filename,clobber=True)
    
    def dumpConfig(self):
        """Dumps a valid configuration file on top of any old configuration files. This is useful for examining the default configuration fo the system, and providing modifications through this interface."""
        with open(self.config["System"]["Configs"]["Source"].rstrip(".yaml")+".dump.yaml",'w') as stream:
            yaml.dump(self.defaults[2]["Source"],stream,default_flow_style=False)


      
