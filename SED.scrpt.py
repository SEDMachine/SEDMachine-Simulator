#!/usr/bin/env python
#
#  SED.scrpt.py
#  Simulation Software
#
#  Created by Alexander Rudy on 2011-11-09.
#  Copyright 2011 Alexander Rudy. All rights reserved.
#

import math, copy, sys, time, logging, os, argparse

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
This is the Command Line Interface to the SEDMachine Simulator Program.

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
        self.ImageDirectory = "Images/"
        self.CacheDirectory = "Caches/"
        self.startLogging()
        self.initOptions()
        self.initSpectra()
        self.bar = arpytools.progressbar.ProgressBar(color="green")
        self.prog = Value('d',0)
    
    def initOptions(self):
        """Set up the options for the command line interface."""
        
        # Establish an argument parser
        self.parser = argparse.ArgumentParser(description=ShortHelp,epilog=LongHelp,
            formatter_class=argparse.RawDescriptionHelpFormatter,)
        
        self.parser.add_argument('title',metavar='label',type=str,help="Include the given string in the filename")
        # Add the basic controls for the script
        self.parser.add_argument('--version',action='version',version=__version__)
        self.parser.add_argument('--dev',action='store_true',dest='dev',help="Enable Development Mode Settings")
        self.parser.add_argument('--no-cache',action='store_false',dest='cache',
            help="Re-run all algorithm parts")
        self.parser.add_argument('--thread',action='store_true',dest='thread',
            help="Enable Threading")
        self.parser.add_argument('-n',action='store',metavar='N',type=int,
            help="Limit to the first N spectra")
        self.parser.add_argument('-o',action='store',metavar='N',type=int,
            help="Offset the starting spectra number by N")
        self.parser.add_argument('--debug',action='store_true',dest='debug',help="Enable Debugging Mode")
        self.parser.add_argument('--test',action='store_true',dest='test',help="Enable Test Mode Settings")
        self.parser.add_argument('--config',action='store',dest='config',type=str,default="SED.config.yaml",
            help="Name of the configuration file for the model")
        self.parser.add_argument('--easy',action='store_true',dest='easy',help="Set easy settings for test runs")
    
    def initSpectra(self):
        """Set up the options for handling single spectra objects"""
        specgroup = self.parser.add_argument_group('Spectra')
        specgroup.add_argument('-s',choices='bfs',help="Select Spectrum: b: BlackBody f:Flat s:File")
        specgroup.add_argument('-T',action='store',dest='Temp',type=float,
            default=5000.0,help="BlackBody Temperature to use")
        specgroup.add_argument('-V',action='store',dest='Value',type=float,
            default=1000.0,help="Flat Spectrum Value")
        specgroup.add_argument('-F',action='store',dest='Filename',type=str,
            default="",help="Spectrum Data Filename")
        specgroup.add_argument('-A',action='store',dest='PreAmp',type=float,
            default=1.0,help="Pre-amplification for Spectrum")
    
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
        
        self.log.info("Simulator has initilaized")
    
    
    def run(self):
        """Runs the simulator"""
        start = time.clock()
        self.parseOptions()
        self.log.info("System Setup")
        self.setup()
        self.log.info("Generating Source")
        self.generateAllLenslets()
        self.placeAllLenslets()
        self.log.info("Writing Image")
        self.saveFile()
        self.log.debug("Total Simulation took %2.3gs for %d lenslets with caching %s" % (time.clock() - start,
            len(self.lenslets),"enabled" if self.options.cache else "disabled"))
    
    
    def parseOptions(self):
        """Interprests the options system"""
        self.options = self.parser.parse_args()
        if self.options.test:
            self.options.easy = True
            self.options.cache = False
        if self.options.dev:
            self.options.debug = True
            self.options.easy = True
        if self.options.easy:
            self.options.s = "b"
            self.options.Temp = 5000
            self.options.o = 2150
            self.options.n = 5
        
        if self.options.debug:
            self.debug = True
            self.console.setLevel(logging.DEBUG)
        else:
            self.logfile.setLevel(logging.INFO)
        
        
        self.optstring = "\n"
        
        for key,value in vars(self.options).iteritems():
            self.optstring += "%(key)15s : %(value)15s" % { 'key':key , 'value':value }
    
    
    def setup(self):
        """Performs all setup options"""
        self.setupSource()
        self.setupModel()
        self.setupLenslets()
    
    def setupSource(self):
        """Sets up the source"""
        if self.options.s == "b":
            self.Spectrum = BlackBodySpectrum(self.options.Temp)
            # self.Spectrum *= self.options.PreAmp
        elif self.options.s == "f":
            self.Spectrum = FlatSpectrum(self.options.Value)
            # self.Spectrum *= self.options.PreAmp
        elif self.options.spec == "s":
            try:
                WL,Flux = np.genfromtxt(self.options.Filename).T
                WL *= 1e-10
                Flux *= 1e18 * 1e6
                self.Spectrum = AS(np.array([WL,Flux]),self.options.Filename)
                # self.Spectrum *= self.options.PreAmp
            except IOError as e:
                self.parser.error("Cannot find Spectrum File: %s" % str(e))
        self.log.debug("Set Spectrum to %s" % self.Spectrum)
    
    def setupModel(self):
        """Sets up the SED Module Model"""
        if self.debug:
            start = time.clock()
        if self.options.cache:
            self.Model = SED.Model(self.options.config,self.CacheDirectory)
        else:
            self.Model = SED.Model(self.options.config)
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
    
    
    def generateLenslet(self,i,spectrum):
        """Generate a single lenslet spectrum"""
        try:
            self.Model.cache_sed_subimage(i,spectrum)
            
            if self.debug:
                self.Model.show()
                plt.savefig("Images/%04d-Subimage.pdf" % i)
                plt.clf()
        
        except SED.SEDLimits:
            msg = "Skipped Lenslet %d, Limits Out of Bounds" % i
            SED.LOG.info(msg)
        else:
            msg = "Cached Lenslet %d" % i
            SED.LOG.info(msg)
        finally:
            self.prog.value += 1.0
            self.bar.render(int((self.prog.value/self.total) * 100),"L:%4d" % i)
    
    def generateAllLenslets(self):
        """Generate all lenslet spectra"""
        self.log.info("Generating Spectra in %d Lenslets" % len(self.lenslets))
        self.bar.render(0,"L:%d" % 0)
        if self.options.thread:
            map_func = self.pool.async_map
        else:
            map_func = map
        map_func(lambda i: self.generateLenslet(i,self.Spectrum),self.lenslets)
    
    
    def placeLenslet(self,i):
        """Place a single lenslet into the model"""
        try:
            self.Model.place_cached_sed(i,"Included Spectrum %d" % i,"Blank")
            if self.debug:
                self.Model.show()
                plt.savefig("Images/%04d-Fullimage.pdf" % i)
                plt.clf()
        except SED.SEDLimits:
            msg = "Encoutered Spectrum outside image boundaries %d" % i
            SED.LOG.info(msg)
        else:
            msg = "Placed Spectrum %d into image" % i
            SED.LOG.info(msg)
        finally:
            self.prog.value += 1
            self.bar.render(int((self.prog.value/self.total) * 100),"L:%4d" % i)
    
    def placeAllLenslets(self):
        """Place all lenslets into the image file"""
        self.log.info("Placing Spectra in %d Lenslets" % len(self.lenslets))
        self.bar.render(0,"L:%d" % 0)
        self.prog.value = 0
        if self.options.thread:
            map_func = self.pool.async_map
        else:
            map_func = map
        map_func(lambda i: self.placeLenslet(i),self.lenslets)
    
    
    def saveFile(self):
        """Saves the file"""
        self.Filename = "%(label)s-%(date)s.%(fmt)s" % { 'label': self.options.title, 'date': time.strftime("%Y-%m-%d"), 'fmt':'fits' }
        self.Fullname = self.ImageDirectory + self.Filename
        self.Model.write(self.Fullname,clobber=True)
        self.log.info("Wrote %s" % self.Fullname)
    


if __name__ == '__main__':
    Sim = Simulator()
    Sim.run()


