#
#  SED.scrpt.py
#  Simulation Software
#
#  Created by Alexander Rudy on 2011-11-09.
#  Copyright 2011 Alexander Rudy. All rights reserved.
#  Version 0.1.2
#

import math, copy, sys, time, logging, os, argparse, yaml, collections

import logging.handlers

import arpytools.progressbar

from multiprocessing import Pool, Value

import numpy as np

import AstroObject.AstroSimulator
import AstroObject.Utilities as AOU

def update(d, u):
    """A deep update command for dictionaries.
    This is because the normal dictionary.update() command does not handle nested dictionaries."""
    if len(u)==0:
        return d
    for k, v in u.iteritems():
        if isinstance(v, collections.Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d

__version__ = AOU.getVersion(__name__)

LongHelp = """
This is the Command Line Interface to the SEDMachine Simulator Program."""

""""
The program works in the following stages:
    1. Instrument Property Setup
    2. Spectrum Generation
    3. Spectra Sub-Image Generation
    4. Full Image Generation (Placing Sub-Images)
    5. Noise Image Generation
    6. Scattered Light Generation
"""
ShortHelp = "Utilities for simulating the SEDMachine"

class Simulator(AstroObject.AstroSimulator.Simulator):
    """A basic SED simulator controller"""
    def __init__(self):
        super(Simulator, self).__init__()
        self.debug = False
        self.lenslets = []
        self.defaults = []
        self.defaultConfig()
        self.bar = arpytools.progressbar.ProgressBar(color="green")
        self.prog = Value('d',0)
        self.initStages()
        
    
    def initStages(self):
        """Set up the options for the command line interface."""
        self.registerStage(self.setupModel,"instinit",description="Set up the Instrument Object")
        self.registerStage(self.setupLenslets,"lensinit",description="Load and verify lenslet objects")
        self.registerMacro("instrument","instinit","lensinit",help="Test the initialization of the instrument model")
        self.registerStage(self.setupNoise,"noisemask",description="Generate Noise Masks")
        self.registerStage(self.setupSource,"sourceinit",description="Set up the source model")
        self.registerMacro("source","sourceinit",help="Test the initialization of the source model")
        self.registerStage(self.generateAllLenslets,"gensubimg",description="Generate Subimages for all Lenslets")
        self.registerStage(self.positionCaches,"cachesubimg",description="Caches Lenslets")
        self.registerStage(self.placeAllLenslets,"placesubimg",description="Places the lenslets into the image")
        self.registerStage(self.cropImage,"crop",description="Crops the image down to size")
        self.registerStage(self.applyNoise,"addnoise")
        self.registerStage(self.saveFile,"save")
        
        
    def defaultConfig(self):
        """Default configuration values from the program"""
        config = {}
        
        config["System"] = {}
        
        # Caching and Logging
        config["System"]["Cache"] = True
        config["System"]["Debug"] = False
        config["System"]["Plot"]  = False
        
        config["System"]["Output"] = {}
        config["System"]["Output"]["Label"] = "Generated"
        config["System"]["Output"]["Format"] = "fits"
        
        config["System"]["CacheFiles"] = {}
        config["System"]["CacheFiles"]["Instrument"] = "SED.instrument"
        config["System"]["CacheFiles"]["Source"] = "SED.source"
        
        # Configuration Files
        config["System"]["Configs"] = {}
        config["System"]["Configs"]["Instrument"] = "SED.instrument.config.yaml"
        config["System"]["Configs"]["Source"] = "SED.source.config.yaml"
        config["System"]["Configs"]["Main"] = "SED.script.config.yaml"
        
        # Directory Configuration
        config["System"]["Dirs"] = {}
        config["System"]["Dirs"]["Logs"] = "Logs/"
        config["System"]["Dirs"]["Partials"] = "Partials/"
        config["System"]["Dirs"]["Caches"] = "Caches/"
        config["System"]["Dirs"]["Images"] = "Images/"
        
        
        # Lenslet Configuration
        config["System"]["Lenslets"] = {}
        
        # Source Configuration
        config["System"]["Source"] = {}
        
        # Logging Configuration
        config["logging"] = {}
        config["logging"]["console"] = {}
        config["logging"]["console"]["enable"] = True
        config["logging"]["console"]["format"] = "...%(message)s"
        config["logging"]["console"]["level"] = logging.INFO
        config["logging"]["file"] = {}
        config["logging"]["file"]["enable"] = True
        config["logging"]["file"]["filename"] = "SED"
        config["logging"]["file"]["format"] = "%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s"
        
        self.defaults += [copy.deepcopy(config)]
        
        self.debug = config["System"]["Debug"]
        
        update(self.config,config)
        self.log.debug("Set up default configutaion")
    
    
    def setupSource(self):
        """Sets up the source"""
        sourceCfg = {}
        if hasattr(self.options,'s'):
            sourceCfg["PreAmp"] = self.options.PreAmp
            if self.options.s == "n":
                self.parser.error("Source Spec Required! See -s")
            elif self.options.s == "b":
                sourceCfg["Type"] = "BlackBody"
                sourceCfg["Temp"] = self.options.temp
            elif self.options.s == "f":
                sourceCfg["Type"] = "Flat"
                sourceCfg["Value"] = self.options.Value
            elif self.options.spec == "s":
                sourceCfg["Type"] = "File"
                sourceCfg["Filename"] = self.options.Filename
        
        
        update(self.config["System"]["Source"],sourceCfg)
        import Source
        
        self.Source = Source.Source(self.config)
        
        self.Source.setup()
                
        self.log.debug("Set Spectrum to %s" % self.Source.Spectrum)
        
        
    
    def setupModel(self):
        """Sets up the SED Module Model"""
        if self.debug:
            start = time.clock()
        
        import Instrument
        self.Limits = Instrument.SEDLimits
        self.Model = Instrument.Instrument(self.config)
        
        self.Model.setup()
        
        if self.debug:
            import matplotlib.pyplot as plt
            self.plt = plt
            end = time.clock()
            dur = end - start
            msg = "Setup took %1.5gs with caches %s." % (dur,"enabled" if self.options.cache else "disabled")
            self.log.debug(msg)
        
        with open("%(dir)s%(fname)s%(fmt)s" % {'dir': self.config["System"]["Dirs"]["Partials"], 'fname': "Instrument-Audit", 'fmt':".dat" },'w') as stream:
            stream.write("State Audit File %s\n" % (time.strftime("%Y-%m-%d-%H:%M:%S")))
        
    
    def setupLenslets(self):
        """Establish the list of lenslets for use in the system"""
        self.lenslets = self.Model.lenslets
        
        if "start" in self.config["System"]["Lenslets"]:
            self.lenslets = self.lenslets[self.config["System"]["Lenslets"]["start"]:]
        if "number" in self.config["System"]["Lenslets"]:
            self.lenslets = self.lenslets[:self.config["System"]["Lenslets"]["number"]]
        
        self.total = len(self.lenslets)
    
    def setupNoise(self):
        """Sets up the noise masks in the model"""
        self.Model.setupNoise()
        
    def positionCaches(self):
        """Caches the positions"""
        self.Model.positionCaching()
    
    def generateLenslet(self,i,spectrum):
        """Generate a single lenslet spectrum"""
        try:
            self.Model.cache_sed_subimage(i,spectrum,write=self.config["System"]["Cache"])
            
            if self.debug:
                self.Model.show()
                self.plt.savefig(self.config["System"]["Dirs"]["Partials"] + "Subimage-%04d-Final.pdf" % i)
                self.plt.clf()
        
        except self.Limits:
            msg = "Skipped Lenslet %4d, Limits Out of Bounds" % i
            self.Model.log.info(msg)
        else:
            msg = "Cached Lenslet %4d" % i
            self.Model.log.info(msg)
        finally:
            self.prog.value += 1.0
            States = self.Model.list()
            with open("%(dir)s%(fname)s%(fmt)s" % {'dir': self.config["System"]["Dirs"]["Partials"], 'fname': "Instrument-Audit", 'fmt':".dat" },'a') as stream:
                stream.write("%4d: %s\n" % (len(States),States))
            self.log.debug("Memory Status: %d states saved in object" % len(States))
            
        if self.debug:
            self.log.info("=>Finished Generating Spectrum %d" % self.prog.value)
        else:
            self.bar.render(int((self.prog.value/self.total) * 100),"L:%4d" % i)
        
    
    
    def generateAllLenslets(self):
        """Generate all lenslet spectra"""
        self.log.info("Generating Spectra in %d Lenslets" % len(self.lenslets))
        if not self.debug:
            self.bar.render(0,"L:%4d" % 0)
        map(lambda i: self.generateLenslet(i,self.Source.getSpectrum(i)),self.lenslets)
        self.bar.lines = 0
    
    
    def placeLenslet(self,i):
        """Place a single lenslet into the model"""
        try:
            self.Model.place_cached_sed(i,"Included Spectrum %d" % i,"Final",fromfile=self.config["System"]["Cache"])
            if self.debug:
                self.Model.show()
                self.plt.savefig(self.config["System"]["Dirs"]["Partials"] + "FullImage-%04d-Final.pdf" % i)
                self.plt.clf()
        except self.Limits:
            msg = "Encoutered Spectrum outside image boundaries %d" % i
            self.Model.log.info(msg)
        else:
            msg = "Placed Spectrum %d into image" % i
            self.Model.log.info(msg)
        finally:
            self.prog.value += 1
            States = self.Model.list()
            with open("%(dir)s%(fname)s%(fmt)s" % {'dir': self.config["System"]["Dirs"]["Partials"], 'fname': "Instrument-Audit", 'fmt':".dat" },'a') as stream:
                stream.write("%4d: %s\n" % (len(States),States))
            self.log.debug("Memory Status: %d states saved in object" % len(States))
            
        if self.debug:
            self.log.info("=>Finished Placing Spectrum %d" % self.prog.value)
        else:
            self.bar.render(int((self.prog.value/self.total) * 100),"L:%4d" % i)
    
    def placeAllLenslets(self):
        """Place all lenslets into the image file"""
        self.Model.save(self.Model.data("Blank"),"Final")
        self.log.info("Placing Spectra in %d Lenslets" % len(self.lenslets))
        if not self.debug:
            self.bar.render(0,"L:%4d" % 0)
        self.prog.value = 0
        map(lambda i: self.placeLenslet(i),self.lenslets)
        self.bar.lines = 0
    
    def cropImage(self):
        """Crops the image down to size"""
        self.Model.ccdCrop()
    
    def applyNoise(self):
        """Applies the noise masks to the file"""
        self.Model.applyNoise(self.Model.statename)
    
    def saveFile(self):
        """Saves the file"""
        self.Filename = "%(label)s-%(date)s.%(fmt)s" % { 'label': self.config["System"]["Output"]["Label"], 'date': time.strftime("%Y-%m-%d"), 'fmt':self.config["System"]["Output"]["Format"] }
        self.Fullname = self.config["System"]["Dirs"]["Images"] + self.Filename
        self.Model.keep(self.Model.statename)
        self.Model.write(self.Fullname,clobber=True)
        self.log.info("Wrote %s" % self.Fullname)
    
def run():
    """Run this simulator"""
    SIM = Simulator()
    SIM.run()
