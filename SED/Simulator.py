#
#  SED.scrpt.py
#  Simulation Software
#
#  Created by Alexander Rudy on 2011-11-09.
#  Copyright 2011 Alexander Rudy. All rights reserved.
#  Version 0.3.0
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
        self.registerStage(self.setupModel,"instinit",help="Set up instrument object",description="Set up the instrument object")
        self.registerMacro("instrument","instinit",help="Test the initialization of the instrument model")
        self.registerStage(self.debugLenslets,"lensdebug",description="Debugging Lenslets",include=False)
        self.registerStage(self.setupNoise,"noisemask",help="Generate noise frames",description="Generate noise masks")
        self.registerStage(self.setupSource,"sourceinit",description="Set up the source model")
        self.registerStage(self.flatSource,"flatsource",description="Set up the flat source model")
        self.registerMacro("source","sourceinit",help="Test the initialization of the source model")
        self.registerStage(self.dispersion,"dispers",description="Get Dispersion for all lenslets")
        self.registerMacro("dispersion","instrument","dispers",help="Generate dispersions for all lenslets")
        self.registerStage(self.trace,"tracer",description="Get Trace for all lenslets")
        self.registerMacro("trace","instrument","dispers","tracer","source",help="Generate dispersions for all lenslets")
        self.registerStage(self.lenslet_image,"subimg",description="Generate subimages for all lenslets")
        self.registerMacro("subimages","instinit","lensinit","sourceinit","dispers","tracer","subimg",help="Generate Sub Images")
        self.registerStage(self.positionCaches,"cachesubimg",description="Cache lenslets")
        self.registerStage(self.placeAllLenslets,"placesubimg",description="Place the lenslets into the image")
        self.registerStage(self.cropImage,"crop",description="Crop the image down to size")
        self.registerStage(self.applyNoise,"addnoise",description="Add noise frames")
        self.registerStage(self.saveFile,"save",description="Save the final image")
        
        self.Caches.registerCustom("Config",kind=AstroObject.AstroSimulator.YAMLCache,filename="Simulator.cache.yaml",generate=lambda : self.config)
        
        self.registerStage(self.plotNoise,"noiseplot",help="3D Plot of Noise Levels",description="3D Plot of Noise Levels",include=False)
        
    def defaultConfig(self):
        """Default configuration values from the program"""
        config = {}
                
        # Caching and Logging
        config["Cache"] = True
        config["Debug"] = False
        config["Plot"]  = False
        
        config["Output"] = {}
        config["Output"]["Label"] = "Generated"
        config["Output"]["Format"] = "fits"
        
        config["CacheFiles"] = {}
        config["CacheFiles"]["Instrument"] = "SED.instrument"
        config["CacheFiles"]["Source"] = "SED.source"
        
        # Configuration Files
        config["Configs"] = {}
        config["Configs"]["Instrument"] = "SED.instrument.config.yaml"
        config["Configs"]["Source"] = "SED.source.config.yaml"
        config["Configs"]["Main"] = "SED.script.config.yaml"
        config["Configs"]["This"] = config["Configs"]["Main"]
        
        # Directory Configuration
        config["Dirs"] = {}
        config["Dirs"]["Logs"] = "Logs/"
        config["Dirs"]["Partials"] = "Partials/"
        config["Dirs"]["Caches"] = "Caches/"
        config["Dirs"]["Images"] = "Images/"
        
        
        # Lenslet Configuration
        config["Lenslets"] = {}
        
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
        
        self.debug = config["Debug"]
        
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
        
        
        update(self.config,sourceCfg)
        import Source
        
        self.Source = Source.Source(self.config)
        
        self.Source.setup()
                
        self.log.debug("Set Spectrum to %s" % self.Source.Spectrum)
        
    def flatSource(self):
        """Make a flat source"""
        self.Source.Spectrum = self.Source.D_Spectrum
    
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
        
    def debugLenslets(self):
        self.Model.plot_lenslet_data()
    
    def setupNoise(self):
        """Sets up the noise masks in the model"""
        self.Model.setupNoise()
        
    def positionCaches(self):
        """Caches the positions"""
        self.Model.positionCaching()
    
    def dispersion(self):
        """Gets the dispersion for each lenslet"""
        self.Model.lenslet_dispersion()

    def trace(self):
        """Calculate the trace for all lenslets"""
        self.Model.lenslet_trace(self.Source.getSpectrum)
        
    
    def generateLenslet(self,i,spectrum):
        """Generate a single lenslet spectrum"""
        try:
            self.Model.cache_sed_subimage(i,spectrum,write=self.config["Cache"])
            
            if self.debug:
                self.Model.show()
                self.plt.savefig(self.config["Dirs"]["Partials"] + "Subimage-%04d-Final.pdf" % i)
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
            with open("%(dir)s%(fname)s%(fmt)s" % {'dir': self.config["Dirs"]["Partials"], 'fname': "Instrument-Audit", 'fmt':".dat" },'a') as stream:
                stream.write("%4d: %s\n" % (len(States),States))
            self.log.debug("Memory Status: %d states saved in object" % len(States))
            
        if self.debug:
            self.log.info("=>Finished Generating Spectrum %d" % self.prog.value)
        else:
            self.bar.render(int((self.prog.value/self.total) * 100),"L:%4d" % i)
        
    
    
    def generateAllLenslets(self):
        """Generate all lenslet spectra"""
        self.Model.cache_images(self.Source.getSpectrum)    
    
    def lenslet_image(self):
        """docstring for generateSubimages"""
        self.Model.lenslet_image()
    
    def placeLenslet(self,i):
        """Place a single lenslet into the model"""
        try:
            self.Model.place_cached_sed(i,"Included Spectrum %d" % i,"Final",fromfile=self.config["Cache"])
            if self.debug:
                self.Model.show()
                self.plt.savefig(self.config["Dirs"]["Partials"] + "FullImage-%04d-Final.pdf" % i)
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
            with open("%(dir)s%(fname)s%(fmt)s" % {'dir': self.config["Dirs"]["Partials"], 'fname': "Instrument-Audit", 'fmt':".dat" },'a') as stream:
                stream.write("%4d: %s\n" % (len(States),States))
            self.log.debug("Memory Status: %d states saved in object" % len(States))
            
        if self.debug:
            self.log.info("=>Finished Placing Spectrum %d" % self.prog.value)
        else:
            self.bar.render(int((self.prog.value/len(self.Model.lenslets)) * 100),"L:%4d" % i)
    
    def placeAllLenslets(self):
        """Place all lenslets into the image file"""
        self.Model.save(self.Model.data("Blank"),"Final")
        self.log.info("Placing Spectra in %d Lenslets" % len(self.Model.lenslets))
        self.bar = arpytools.progressbar.ProgressBar(color="Red")
        if not self.debug:
            self.bar.render(0,"L:%4d" % 0)
        self.prog.value = 0
        self.Model.log.useConsole(False)
        self.log.useConsole(False)
        map(lambda i: self.placeLenslet(i),self.Model.lenslets)
        self.log.useConsole(True)
        self.Model.log.useConsole(True)
        
        
        self.bar.lines = 0
    
    def cropImage(self):
        """Crops the image down to size"""
        self.Model.ccdCrop()
    
    def applyNoise(self):
        """Applies the noise masks to the file"""
        self.Model.applyNoise(self.Model.statename)
    
    def saveFile(self):
        """Saves the file"""
        self.Filename = "%(label)s-%(date)s.%(fmt)s" % { 'label': self.config["Output"]["Label"], 'date': time.strftime("%Y-%m-%d"), 'fmt':self.config["Output"]["Format"] }
        self.Fullname = self.config["Dirs"]["Images"] + self.Filename
        self.Model.keep(self.Model.statename)
        self.Model.write(self.Fullname,clobber=True)
        self.log.info("Wrote %s" % self.Fullname)
        
    def plotNoise(self):
        """Plot noise"""
        self.Model.show3D("Dark")
        format = { 'dir': self.config["Dirs"]["Partials"], 'stage':"Dark", 'fmt':self.config["plot_format"]}
        filename = "%(dir)s/%(stage)s.%(fmt)s" 
        plt.savefig(filename % format)
        plt.clf()
        self.Model.show3D("Bias")
        format["stage"] = "Bias"
        plt.savefig(filename % format)
        plt.clf()
    
def run():
    """Run this simulator"""
    SIM = Simulator()
    SIM.run()
