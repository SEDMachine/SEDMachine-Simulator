#!/usr/bin/env python
#
#  SED.scrpt.py
#  Simulation Software
#
#  Created by Alexander Rudy on 2011-11-09.
#  Copyright 2011 Alexander Rudy. All rights reserved.
#

import math, copy, sys, time, logging, os, argparse, yaml

import numpy as np
import pyfits as pf
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

from multiprocessing import Pool, Value

import SED
from AstroObject.AnalyticSpectra import BlackBodySpectrum,FlatSpectrum,InterpolatedSpectrum
from AstroObject.AstroSpectra import SpectraObject
import arpytools.progressbar
import scipy.signal

__version__ = SED.__version__

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

class Simulator(object):
    """A basic SED simulator controller"""
    def __init__(self):
        super(Simulator, self).__init__()
        self.debug = False
        self.startLogging()
        self.initOptions()
        self.bar = arpytools.progressbar.ProgressBar(color="green")
        self.prog = Value('d',0)
    
    def initOptions(self):
        """Set up the options for the command line interface."""
        
        USAGE = """SEDrun.py [-D | -T | -E | -F] [arguments] subcommand"""
        # Establish an argument parser
        self.parser = argparse.ArgumentParser(description=ShortHelp,epilog=LongHelp,
            formatter_class=argparse.RawDescriptionHelpFormatter,usage=USAGE)
        
        self.parser.add_argument('-f',metavar='label',type=str,dest='title',
            help="label for output image",default="Experiment")
        # Add the basic controls for the script
        self.parser.add_argument('--version',action='version',version=__version__)
        
        self.basics = self.parser.add_mutually_exclusive_group()
        self.basics.add_argument('-D','--dev',action='store_true',dest='dev',
            help="equivalent to --debug --easy")
        self.basics.add_argument('-T','--test',action='store_true',dest='test',
            help="equivalent to --no-cache --easy")
        self.basics.add_argument('-E','--easy',action='store_true',dest='easy',
            help="uses simple configuration for testing the simulator")
        self.basics.add_argument('-F','--full',action='store_true',dest='full',
            help="disables --debug and forces normal deployment operation")
            
        self.parser.add_argument('--plot',action='store_true',dest='plot',help="Enable debugging plots")
        self.parser.add_argument('--no-cache',action='store_false',dest='cache',
            help="ignore cached data from the simulator")
        self.parser.add_argument('--thread',action='store_true',dest='thread',
            help="enable multi-process support")
        self.parser.add_argument('-n',action='store',metavar='N',type=int,
            help="limit the run to N lenslets")
        self.parser.add_argument('-o',action='store',metavar='I',type=int,
            help="start filling from lenslet in array position I")
            
        self.parser.add_argument('-d','--debug',action='store_true',dest='debug',help="enable debugging messages and plots")
        
        self.parser.add_argument('--config',action='store',dest='config',type=str,default="SED.script.config.yaml",
            help="use the specified configuration file",metavar="file.yaml")
        self.parser.add_argument('--iconfig',action='store',dest='iconfig',type=str,default="SED.instrument.config.yaml",
            help="use the specified configuration file",metavar="file.yaml")
        self.parser.add_argument('--dump-config',action='store_true',dest='dump',help=argparse.SUPPRESS)
        
        
        self.uniform = argparse.ArgumentParser(add_help=False)
        
        self.uniform.add_argument('-s',choices='bfs',default='n',help="Select Spectrum: b: BlackBody f:Flat s:File")
        self.uniform.add_argument('-T',action='store',dest='Temp',type=float,
            default=5000.0,help="BlackBody Temperature to use")
        self.uniform.add_argument('-V',action='store',dest='Value',type=float,
            default=1000.0,help="Flat Spectrum Value")
        self.uniform.add_argument('-F',action='store',dest='Filename',type=str,
            default="",help="Spectrum Data Filename")
        self.uniform.add_argument('-A',action='store',dest='PreAmp',type=float,
            default=1.0,help="Pre-amplification for Spectrum")
        
        self.subparsers = self.parser.add_subparsers(title="sub-commands",dest="command",metavar="subcommand")
        
        self.initSpectra()
        self.initStartup()
        self.initSource()
        self.initPosTest()
    
    def initSpectra(self):
        """Set up the options for handling single spectra objects"""
        ShortHelp = "create images with a uniform source field spectrum"
        Description = "Runs the simulator with a uniform source, i.e. the same spectrum in each lenslet"
        specgroup = self.subparsers.add_parser('uniform',description=Description,help=ShortHelp,parents=[self.uniform])
    
    def initPosTest(self):
        """Position Testing Script"""
        ShortHelp = "position testing mode for spectrum placement"
        Description = "Produces output to determine if the spectrum placement appears to have been done correctly"
        postest = self.subparsers.add_parser('postest',description=Description,help=ShortHelp,parents=[self.uniform])
    
    def initStartup(self):
        """Subcommand for handling only startup functions"""
        ShortHelp = "run the initalization and caching for the system"
        Description = "Initializes the simulation."
        startupgroup = self.subparsers.add_parser('startup',description=Description,help=ShortHelp)
    
    def initSource(self):
        """Subcommand for handling only startup and source functions"""
        ShortHelp = "run the source creation and resolution routines"
        Description = "Initializes the simulation, then initializes the source spectra"
        sourcegroup = self.subparsers.add_parser('source',description=Description,help=ShortHelp,parents=[self.uniform])
    
    def startLogging(self):
        """Establishes logging for this module"""
        
        self.log = logging.getLogger(sys.argv[0])
        logfolder = "Logs/"
        filename = sys.argv[0]+"-"+time.strftime("%Y-%m-%d")
        longFormat = "%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s"
        shortFormat = '... %(message)s'
        dateFormat = "%Y-%m-%d-%H:%M:%S"
        
        self.log.setLevel(logging.DEBUG)
        logging.captureWarnings(True)
        self.console = logging.StreamHandler()
        self.console.setLevel(logging.INFO)
        consoleFormatter = logging.Formatter(shortFormat,datefmt=dateFormat)
        self.console.setFormatter(consoleFormatter)
        self.log.addHandler(self.console)
        
        self.logfile = logging.FileHandler(filename=logfolder+filename+".log",mode="a")
        self.logfile.setLevel(logging.DEBUG)
        fileformatter = logging.Formatter(longFormat,datefmt=dateFormat)
        self.logfile.setFormatter(fileformatter)
        
        if os.access(logfolder,os.F_OK):
            self.log.addHandler(self.logfile)
        
        self.log.info("Runner has initilaized")
    
    
    def run(self):
        """Runs the simulator"""
        start = time.clock()
 
        self.parseOptions()
        cmd = self.options.command
        
        if self.options.dump:
            self.log.info("Dumping Configuration Variables...")
            self.dumpConfig()
            self.log.info("Cancelling all subcommands")
            cmd = ""
        
        if cmd in ["uniform","source","startup","postest"]:
            self.log.info("Simulator Setup")
            self.setup()
        
        if cmd in ["uniform","source","postest"]:
            self.log.info("Source Setup")
            self.setupSource()
        
        if cmd in ["postest"]:
            self.log.info("Testing Source Positioning")
            self.positionTests()
        
        if cmd in ["uniform"]:
            self.log.info("Generating Source")
            self.generateAllLenslets()
            
        if cmd in ["uniform"]:
            self.placeAllLenslets()
            self.log.info("Writing Image")
            
        if cmd in ["uniform"]:
            self.saveFile()
        
        if cmd in ["uniform","source","postest"]:
            self.log.debug("Total Simulation took %2.3gs for %d lenslets with caching %s" % (time.clock() - start,
                len(self.lenslets),"enabled" if self.options.cache else "disabled"))
        self.log.info("Runner is done")
    
    
    def parseOptions(self):
        """Interprests the options system"""
        self.options = self.parser.parse_args()
        
        msg = ""
        
        if self.options.test:
            self.options.easy = True
            self.options.cache = False
            msg = "Test Mode"
        if self.options.dev:
            self.options.debug = True
            self.options.easy = True
            msg = "Dev Mode"
        if self.options.easy:
            msg += "Easy Settings"
            if not hasattr(self.options,'s') or self.options.s == "n":
                self.options.s = "b"
                self.options.Temp = 5000
            if not self.options.o:
                self.options.o = 2150
            if not self.options.n:
                self.options.n = 5
        
        if msg != "":
            self.log.info("Registering %s" % msg)
        
        if self.options.debug:
            self.debug = True
            self.console.setLevel(logging.DEBUG)
        
        if self.options.thread:
            self.parser.error("Threading is not yet implemented: AstroObjects are not Thread-Safe")
        
        self.loadConfig()
        
        self.optstring = "SEDScript: \n"
        
        for key,value in vars(self.options).iteritems():
            self.optstring += "%(key)15s : %(value)15s \n" % { 'key':key , 'value':value }
        
        stream = file(self.config["Dirs"]["Partials"]+"/Options.dat",'w')
        stream.write(self.optstring)
        stream.close()
    
    def loadConfig(self):
        """Loads the config file"""
        self.config = {"Dirs":{"Caches":"Caches/","Images":"Images/","Partials":"Partials/","Logs":"Logs/"},"Cache":self.options.cache}
        try:
            self.config = SED.update(self.config,yaml.load(file(self.options.config)))
        except IOError:
            self.log.warning("Couldn't Read configuration file %s, using defaults" % self.options.config)
        self.dirCheck()
        
    def dirCheck(self):
        """Checks that all of the configured directories exist, and creates the ones that don't."""
        for DIR in self.config["Dirs"].values():
            if not os.access(DIR,os.F_OK):
                os.mkdir(DIR)
                self.log.info("Created Directory %s" % DIR)
        
    def dumpConfig(self):
        """Dumps a config back out"""
        yaml.dump(self.config,file(self.options.config,'w'),default_flow_style=False)
        Model = SED.Model(self.options.iconfig,self.config)
        Model.dumpConfig()
        
    def setup(self):
        """Performs all setup options"""
        self.setupModel()
        self.setupLenslets()
    
    def setupSource(self):
        """Sets up the source"""
        if self.options.s == "n":
            self.parser.error("Source Spec Required! See -s")
        elif self.options.s == "b":
            self.Spectrum = BlackBodySpectrum(self.options.Temp)
            self.Spectrum *= self.options.PreAmp
        elif self.options.s == "f":
            self.Spectrum = FlatSpectrum(self.options.Value)
            self.Spectrum *= self.options.PreAmp
        elif self.options.spec == "s":
            try:
                WL,Flux = np.genfromtxt(self.options.Filename).T
            except IOError as e:
                self.parser.error("Cannot find Spectrum File: %s" % str(e))
            WL *= 1e-10
            Flux *= 1e18 * 1e6
            self.Spectrum = AS(np.array([WL,Flux]),self.options.Filename)
            self.Spectrum *= self.options.PreAmp
        self.log.debug("Set Spectrum to %s" % self.Spectrum)
        
        
    
    def setupModel(self):
        """Sets up the SED Module Model"""
        if self.debug:
            start = time.clock()
        
        self.Model = SED.Model(self.options.iconfig,self.config)
        
        self.Model.plot = self.options.plot
        self.Model.setup()
        
        if self.debug:
            end = time.clock()
            dur = end - start
            msg = "Setup took %1.5gs with caches %s." % (dur,"enabled" if self.options.cache else "disabled")
            self.log.debug(msg)
    
    def setupLenslets(self):
        """Establish the list of lenslets for use in the system"""
        self.lenslets = self.Model.lenslets
        
        if self.options.o:
            self.lenslets = self.lenslets[self.options.o:]
        if self.options.n:
            self.lenslets = self.lenslets[:self.options.n]
        
        self.total = len(self.lenslets)
    
    
    def positionTests(self):
        """Test the positioning of spectra on the image"""
        self.bar.render(0,"L:%4d" % 0)
        handle = file(self.config["Dirs"]["Partials"] + "Positions.dat",'w')
        handle.write("# Spectra Positions\n")
        handle.close()
        for i in self.lenslets:
            self.positionTest(i,self.Spectrum)
        self.bar.lines = 0
    
    def positionTest(self,lenslet,spectrum):
        """A single position test"""
        try:
            image,corner = self.Model.get_dense_image(lenslet,spectrum)
            points,wl,deltawl = self.Model.get_wavelengths(lenslet)
            x,y = points.T
            ncorner = np.array([np.max(x),np.min(y)])
            handle = file(self.config["Dirs"]["Partials"] + "Positions.dat",'a')
            np.savetxt(handle,np.array([np.hstack((corner,ncorner))]),fmt='%6.1f')
        except SED.SEDLimits:
            msg = "Skipped Lenslet %4d" % lenslet
            self.Model.log.info(msg)
        else:
            msg = "Cached Lenslet %4d" % lenslet
            self.Model.log.info(msg)
        finally:
            self.prog.value += 1.0
            self.bar.render(int((self.prog.value/self.total) * 100),"L:%4d" % lenslet)
    
    def generateLenslet(self,i,spectrum):
        """Generate a single lenslet spectrum"""
        try:
            self.Model.cache_sed_subimage(i,spectrum)
            
            if self.debug:
                self.Model.show()
                plt.savefig(self.config["Dirs"]["Partials"] + "%04d-Subimage.pdf" % i)
                plt.clf()
        
        except SED.SEDLimits:
            msg = "Skipped Lenslet %4d, Limits Out of Bounds" % i
            self.Model.log.info(msg)
        else:
            msg = "Cached Lenslet %4d" % i
            self.Model.log.info(msg)
        finally:
            self.prog.value += 1.0
            self.bar.render(int((self.prog.value/self.total) * 100),"L:%4d" % i)
    
    def generateAllLenslets(self):
        """Generate all lenslet spectra"""
        self.log.info("Generating Spectra in %d Lenslets" % len(self.lenslets))
        self.bar.render(0,"L:%4d" % 0)
        if self.options.thread:
            map_func = self.pool.async_map
        else:
            map_func = map
        map_func(lambda i: self.generateLenslet(i,self.Spectrum),self.lenslets)
        self.bar.lines = 0
    
    
    def placeLenslet(self,i):
        """Place a single lenslet into the model"""
        try:
            self.Model.place_cached_sed(i,"Included Spectrum %d" % i,"Final")
            if self.debug:
                self.Model.show()
                plt.savefig(self.config["Dirs"]["Partials"] + "%04d-Fullimage.pdf" % i)
                plt.clf()
        except SED.SEDLimits:
            msg = "Encoutered Spectrum outside image boundaries %d" % i
            self.Model.log.info(msg)
        else:
            msg = "Placed Spectrum %d into image" % i
            self.Model.log.info(msg)
        finally:
            self.prog.value += 1
            self.bar.render(int((self.prog.value/self.total) * 100),"L:%4d" % i)
    
    def placeAllLenslets(self):
        """Place all lenslets into the image file"""
        self.Model.save(self.Model.data("Blank"),"Final")
        self.log.info("Placing Spectra in %d Lenslets" % len(self.lenslets))
        self.bar.render(0,"L:%4d" % 0)
        self.prog.value = 0
        if self.options.thread:
            map_func = self.pool.async_map
        else:
            map_func = map
        map_func(lambda i: self.placeLenslet(i),self.lenslets)
        self.bar.lines = 0
    
    def saveFile(self):
        """Saves the file"""
        self.Filename = "%(label)s-%(date)s.%(fmt)s" % { 'label': self.options.title, 'date': time.strftime("%Y-%m-%d"), 'fmt':'fits' }
        self.Fullname = self.config["Dirs"]["Images"] + self.Filename
        self.Model.keep(self.Model.statename)
        self.Model.write(self.Fullname,clobber=True)
        self.log.info("Wrote %s" % self.Fullname)
    


if __name__ == '__main__':
    Sim = Simulator()
    Sim.run()


