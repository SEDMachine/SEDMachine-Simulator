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
        self.registerStage(self.setupModel,"Instrument Setup",description="Set up the Instrument Object")
        self.registerStage(self.setupLenslets,"Lenslet Initalization",description="Load and verify lenslet objects")
        self.registerMacro("instrument","Instrument Setup","Lenslet Initalization",help="Test the initialization of the instrument model")
        self.registerStage(self.setupNoise,"Noise Masks",description="Generate Noise Masks")
        self.registerStage(self.setupSource,"Setup Source",description="Set up the source model")
        self.registerMacro("source","Setup Source")
        self.registerStage(self.generateAllLenslets,"Generate Lenslets",description="Generate Subimages for all Lenslets")
        
    
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
        config["System"]["logging"] = {}
        config["System"]["logging"]["console"] = {}
        config["System"]["logging"]["console"]["enable"] = True
        config["System"]["logging"]["console"]["format"] = "...%(message)s"
        config["System"]["logging"]["console"]["level"] = False
        config["System"]["logging"]["file"] = {}
        config["System"]["logging"]["file"]["enable"] = True
        config["System"]["logging"]["file"]["filename"] = "SED"
        config["System"]["logging"]["file"]["format"] = "%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s"
        
        self.defaults += [copy.deepcopy(config)]
        
        self.debug = config["System"]["Debug"]
        
        update(self.config,config)
        self.log.debug("Set up default configutaion")
    
    def parseOptions(self):
        """Interprests the options system"""
        self.options = self.parser.parse_args()

        self.config["System"]["Configs"]["This"] = self.options.config

        mode = self.options.mode
        
        if mode == None:
            mode = "full"
        
        self.log.debug("Operating Mode: %s" % mode)
        
        if mode != "full":
            self.config["System"]["Lenslets"]["start"] = 2150
            self.config["System"]["Lenslets"]["number"] = 5
            self.config["System"]["Source"]["Type"] = "BlackBody"
            self.config["System"]["Source"]["Temp"] = 5000
        if mode == "dev":
            self.config["System"]["Debug"] = True
            self.config["System"]["Cache"] = True
        elif mode == "test":
            self.config["System"]["Debug"] = False
            self.config["System"]["Cache"] = False
        elif mode == "easy":
            self.config["System"]["Debug"] = False
            self.config["System"]["Cache"] = True

        self.loadConfig(self.config["System"]["Configs"]["This"],self.config["System"])
        
        if hasattr(self.options,'sourceconfig') and type(self.options.sourceconfig)==str:
            self.log.debug("Setting source config to %s" % self.options.sourceconfig)
            self.config["System"]["Configs"]["Source"] = self.options.sourceconfig
        
        if self.options.n:
            self.config["System"]["Lenslets"]["number"] = self.options.n
        
        if self.options.o:
            self.config["System"]["Lenslets"]["start"] = self.options.o
        
        # The trigger sets this to False, as such the option should unset it.
        self.config["System"]["Cache"] &= self.options.cache
        
        # Default value for this is False, as such the option should trigger it.
        self.config["System"]["Debug"] |= self.options.debug
        
        # The option for this turns it on, as such, this should be or.
        self.config["System"]["Plot"] |= self.options.plot
        
        if hasattr(self.options,'label'):
            self.config["System"]["Output"]["Label"] = self.options.label
        
        self.debug = self.config["System"]["Debug"]
        
        self.setupLog()
        
        if self.config["System"]["Debug"]:
            self.console.setLevel(logging.DEBUG)
        
        self.dirCheck()
        
        self.optstring = "SEDScript: \n"

        for key,value in vars(self.options).iteritems():
            self.optstring += "%(key)15s : %(value)-40s \n" % { 'key':key , 'value':value }

        with open(self.config["System"]["Dirs"]["Partials"]+"/Script-Options.dat",'w') as stream:
            stream.write(self.optstring)
    
    def oldrun(self):
        """Runs the simulator"""
        self.start = time.clock()
        
        try:
            self.parseOptions()
            cmd = self.options.command
        
            if self.options.dump:
                self.log.info("Dumping Configuration Variables...")
                self.dumpConfig()
                self.exit()
                return
            
            if cmd in ["postest"]:
                self.log.info("Testing Source Positioning")
                self.positionTests()
                self.exit()
                return
        
            self.log.info("Caching Wavelength Data")
            self.positionCaches()
        
            if cmd in ["cache"]:
                self.exit()
                return
        
            # Log Message generated by subfunction which understands how many lenslets
            self.placeAllLenslets()
        
            self.log.info("Cropping Image")
            self.cropImage()
        
            self.log.info("Applying Noise Masks")
            self.applyNoise()
        
            self.log.info("Writing Image")
            self.saveFile()
        
            self.exit()
            return
        except Exception as e:
            self.log.critical("Simulator Encoutered a Critical Error, and was forced to close!")
            self.log.critical("Exception %s" % str(e))
            raise
        
    
    def loadConfig(self,filename,dest):
        """Loads the config file"""
        try:
            with open(filename) as stream:
                cfg = update(dest,yaml.load(stream))
        except IOError:
            self.log.warning("Couldn't Read configuration file %s, using defaults" % filename)
            cfg = dest
        else:
            self.log.debug("Loaded configuration from file %s" % filename)
        
        return cfg
        
    def dirCheck(self):
        """Checks that all of the configured directories exist, and creates the ones that don't."""
        for DIR in self.config["System"]["Dirs"].values():
            if not os.access(DIR,os.F_OK):
                os.mkdir(DIR)
                self.log.info("Created Directory %s" % DIR)
        
    def dumpConfig(self):
        """Dumps a config back out"""
        import Instrument, Source
        with open(self.config["System"]["Configs"]["This"].rstrip(".yaml")+".dump.yaml",'w') as stream:
            yaml.dump(self.config["System"],stream,default_flow_style=False)
        
        Model = Instrument.Instrument(self.config)
        Model.setup()
        Model.dumpConfig()
        
        Source = Source.Source(self.config)
        Source.setup()
        Source.dumpConfig()
    
    
    
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
        
    
    def positionTests(self):
        """Test the positioning of spectra on the image"""
        self.bar.render(0,"L:%4d" % 0)
        handle = file(self.config["System"]["Dirs"]["Partials"] + "Instrument-Positions.dat",'w')
        handle.write("# Spectra Positions\n")
        handle.close()
        for i in self.lenslets:
            self.positionTest(i,self.Spectrum)
        self.bar.lines = 0
    
    def positionCaches(self):
        """Caches the positions"""
        self.Model.positionCaching()
    
    def positionTest(self,lenslet,spectrum):
        """A single position test"""
        try:
            image,corner = self.Model.get_dense_image(lenslet,spectrum)
            points,wl,deltawl = self.Model.get_wavelengths(lenslet)
            x,y = points.T
            ncorner = np.array([np.max(x),np.min(y)])
            handle = file(self.config["System"]["Dirs"]["Partials"] + "Instrument-Positions.dat",'a')
            np.savetxt(handle,np.array([np.hstack((corner,ncorner))]),fmt='%6.1f')
        except self.Limits:
            msg = "Skipped Lenslet %4d" % lenslet
            self.Model.log.debug(msg)
        else:
            msg = "Cached Lenslet %4d" % lenslet
            self.Model.log.debug(msg)
        finally:
            self.prog.value += 1.0
            self.bar.render(int((self.prog.value/self.total) * 100),"L:%4d" % lenslet)
    
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
        if self.options.thread:
            map_func = self.pool.async_map
        else:
            map_func = map
        map_func(lambda i: self.placeLenslet(i),self.lenslets)
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
