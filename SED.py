#!/usr/bin/env python
#
#  SED.scrpt.py
#  Simulation Software
#
#  Created by Alexander Rudy on 2011-11-09.
#  Copyright 2011 Alexander Rudy. All rights reserved.
#  Version 0.1.1
#

import math, copy, sys, time, logging, os, argparse, yaml, collections

import logging.handlers

import arpytools.progressbar

from multiprocessing import Pool, Value

import numpy as np

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

__version__ = open(os.path.abspath(os.path.join(os.path.dirname(__file__),"VERSION")),'r').read().rstrip("\n")

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
        self.lenslets = []
        self.initLog()
        self.defaults = []
        self.defaultConfig()
        self.bar = arpytools.progressbar.ProgressBar(color="green")
        self.prog = Value('d',0)
        self.initOptions()
        
    
    def initOptions(self):
        """Set up the options for the command line interface."""
        
        self.USAGE = sys.argv[0] + """ [-D | -T | -E | -F] [arguments] %s"""
        # Establish an argument parser
        self.parser = argparse.ArgumentParser(description=ShortHelp,epilog=LongHelp,
            formatter_class=argparse.RawDescriptionHelpFormatter,usage=self.USAGE % "subcommand")
        
        self.parser.add_argument('-f',metavar='label',type=str,dest='title',
            help="label for output image")
        # Add the basic controls for the script
        self.parser.add_argument('--version',action='version',version=__version__)
        
        # Mode Controls
        self.basics = self.parser.add_mutually_exclusive_group()
        self.basics.add_argument('-D','--dev',action='store_const',dest='mode',const='dev',
            help="equivalent to --debug --easy")
        self.basics.add_argument('-T','--test',action='store_const',dest='mode',const='test',
            help="equivalent to --no-cache --easy")
        self.basics.add_argument('-E','--easy',action='store_const',dest='mode',const='easy',
            help="uses simple configuration for testing the simulator")
        self.basics.add_argument('-F','--full',action='store_const',dest='mode',const='full',
            help="disables --debug and forces normal deployment operation")
        
        # Operational Controls
        self.parser.add_argument('--plot',action='store_true',dest='plot',help="Enable debugging plots")
        self.parser.add_argument('--no-cache',action='store_false',dest='cache',
            help="ignore cached data from the simulator")
        self.parser.add_argument('-d','--debug',action='store_true',dest='debug',help="enable debugging messages and plots")
        
        # Special Configuration items
        self.parser.add_argument('-n',action='store',metavar='N',type=int,
            help="limit the run to N lenslets")
        self.parser.add_argument('-o',action='store',metavar='I',type=int,
            help="start filling from lenslet in array position I")
        
        # Config Commands
        self.parser.add_argument('--config',action='store',dest='config',type=str,
            default=self.config["System"]["Configs"]["This"],help="use the specified configuration file",metavar="file.yaml")
        self.parser.add_argument('--dump-config',action='store_true',dest='dump',help=argparse.SUPPRESS)
        
        # Futrure Options
        self.parser.add_argument('--thread',action='store_true',dest='thread',help=argparse.SUPPRESS)
        
        self.all = argparse.ArgumentParser(add_help=False)
        
        self.all.add_argument('-c',action='store',metavar='config',dest='sourceconfig',type=str,
            help="Source configuration file")
        
        self.subparsers = self.parser.add_subparsers(title="sub-commands",dest="command")
        
        self.initAll()
        self.initCache()
        self.initStartup()
        # self.initPosTest()
        self.initSource()
        self.initInstrument()
        
        self.log.debug("Set up Command Line Control Options")
    
    def initAll(self):
        """Set up the options for handling a basic command"""
        Command = "all"
        Usage = self.USAGE % "%s [-c config ]"% (Command)
        ShortHelp = "run the full simulator"
        Description = "Runs the simulator with source configuration provided"
        specgroup = self.subparsers.add_parser(Command,description=Description,help=ShortHelp,
            parents=[self.all],usage=Usage)
    
    def initCache(self):
        """Position Caching Subcommand"""
        Command = "cache"
        Usage = self.USAGE % "%s" % (Command)
        ShortHelp = "Caching mode for spectrum placement - Caches everything except subimages"
        Description = "Caches results of the model._get_wavelenghts results for faster lookup later"
        postest = self.subparsers.add_parser(Command,description=Description,help=ShortHelp,
            usage=Usage)
        
    
    def initPosTest(self):
        """Position Testing Script"""
        Command = "postest"
        Usage = self.USAGE % "%s" % (Command)
        ShortHelp = "position testing mode for spectrum placement"
        Description = "Produces output to determine if the spectrum placement appears to have been done correctly"
        postest = self.subparsers.add_parser(Command,description=Description,help=ShortHelp,usage=Usage)
    
    def initInstrument(self):
        """Subcommand for handling only startup functions"""
        Command = "instrument"
        Usage = self.USAGE % "%s" % (Command)
        ShortHelp = "run the initalization and caching for the instrument"
        Description = "Initializes the instrument."
        startupgroup = self.subparsers.add_parser(Command,description=Description,help=ShortHelp,usage=Usage)
    
    def initStartup(self):
        """Subcommand for handling only startup functions"""
        Command = "startup"
        Usage = self.USAGE % "%s" % (Command)
        ShortHelp = "run the initalization and caching for the system"
        Description = "Initializes the simulation."
        startupgroup = self.subparsers.add_parser(Command,description=Description,help=ShortHelp,usage=Usage)
    
    def initSource(self):
        """Subcommand for handling only startup and source functions"""
        Command = "source"
        Usage = self.USAGE % "%s" % (Command)
        ShortHelp = "run the source creation and resolution routines"
        Description = "Initializes the simulation, then initializes the source spectra"
        sourcegroup = self.subparsers.add_parser(Command,description=Description,help=ShortHelp,parents=[self.all])
    
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
        
        self.log.debug("----------------------------------------")
        self.log.info("Welcome to the SEDMachine Data Simulator")
        self.log.debug("Version %s" % __version__)
        
    def setupLog(self):
        """Setup Logging Functions for the SEDMachine Model.
        
        This configures the logging system, including a possible console and file log. It then reads the logging buffer into the logfile.
        """
        
        # Setup the Console Log Handler
        self.console = logging.StreamHandler()
        consoleFormat = self.config["System"]["logging"]["console"]["format"]
        if self.config["System"]["logging"]["console"]["level"]:
            self.console.setLevel(self.config["System"]["logging"]["console"]["level"])
        elif self.debug:
            self.console.setLevel(logging.DEBUG)
        else:
            self.console.setLevel(logging.INFO)
        consoleFormatter = logging.Formatter(consoleFormat)
        self.console.setFormatter(consoleFormatter)
        
        if self.config["System"]["logging"]["console"]["enable"]:
            self.log.addHandler(self.console)
            self.consolebuffer.setTarget(self.console)
        self.consolebuffer.close()
        self.log.removeHandler(self.consolebuffer)
        
        self.logfile = None
        # Only set up the file log handler if we can actually access the folder
        if os.access(self.config["System"]["Dirs"]["Logs"],os.F_OK) and self.config["System"]["logging"]["file"]["enable"]:
            filename = self.config["System"]["Dirs"]["Logs"] + self.config["System"]["logging"]["file"]["filename"]+".log"
            self.logfile = logging.handlers.TimedRotatingFileHandler(filename=filename,when='midnight')
            self.logfile.setLevel(logging.DEBUG)
            fileformatter = logging.Formatter(self.config["System"]["logging"]["file"]["format"],datefmt="%Y-%m-%d-%H:%M:%S")
            self.logfile.setFormatter(fileformatter)
            self.log.addHandler(self.logfile)
            # Finally, we should flush the old buffers
            self.filebuffer.setTarget(self.logfile)
        
        self.filebuffer.close()
        self.log.removeHandler(self.filebuffer)
        self.log.debug("Configured Logging")
        
    
    def defaultConfig(self):
        """Default configuration values from the program"""
        self.config = {}
        
        self.config["System"] = {}
        
        # Caching and Logging
        self.config["System"]["Cache"] = True
        self.config["System"]["Debug"] = False
        self.config["System"]["Plot"]  = False
        
        self.config["System"]["Output"] = {}
        self.config["System"]["Output"]["Label"] = "Generated"
        self.config["System"]["Output"]["Format"] = "fits"
        
        self.config["System"]["CacheFiles"] = {}
        self.config["System"]["CacheFiles"]["Instrument"] = "SED.instrument"
        self.config["System"]["CacheFiles"]["Source"] = "SED.source"
        
        # Configuration Files
        self.config["System"]["Configs"] = {}
        self.config["System"]["Configs"]["Instrument"] = "SED.instrument.config.yaml"
        self.config["System"]["Configs"]["Source"] = "SED.source.config.yaml"
        self.config["System"]["Configs"]["This"] = "SED.script.config.yaml"
        
        # Directory Configuration
        self.config["System"]["Dirs"] = {}
        self.config["System"]["Dirs"]["Logs"] = "Logs/"
        self.config["System"]["Dirs"]["Partials"] = "Partials/"
        self.config["System"]["Dirs"]["Caches"] = "Caches/"
        self.config["System"]["Dirs"]["Images"] = "Images/"
        
        
        # Lenslet Configuration
        self.config["System"]["Lenslets"] = {}
        
        # Source Configuration
        self.config["System"]["Source"] = {}
        
        # Logging Configuration
        self.config["System"]["logging"] = {}
        self.config["System"]["logging"]["console"] = {}
        self.config["System"]["logging"]["console"]["enable"] = True
        self.config["System"]["logging"]["console"]["format"] = "...%(message)s"
        self.config["System"]["logging"]["console"]["level"] = False
        self.config["System"]["logging"]["file"] = {}
        self.config["System"]["logging"]["file"]["enable"] = True
        self.config["System"]["logging"]["file"]["filename"] = "SED"
        self.config["System"]["logging"]["file"]["format"] = "%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s"
        
        self.defaults += [copy.deepcopy(self.config)]
        
        self.debug = self.config["System"]["Debug"]
        
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

        with open(self.config["System"]["Dirs"]["Partials"]+"/Options.dat",'w') as stream:
            stream.write(self.optstring)
    
    def run(self):
        """Runs the simulator"""
        self.start = time.clock()
        
        self.parseOptions()
        cmd = self.options.command
        
        if self.options.dump:
            self.log.info("Dumping Configuration Variables...")
            self.dumpConfig()
            self.exit()
            return
        
        
        if cmd not in ["source"]:
            self.log.info("Simulator Setup")
            self.setup()
        
        if cmd in ["instrument"]:
            self.exit()
            return
        
        if cmd not in ["source"]:
            self.log.info("Generating Noise Masks")
            self.setupNoise()
        
        self.log.info("Source Setup")
        self.setupSource()
        
        if cmd in ["startup"]:
            self.exit()
            return
        
        if cmd in ["source"]:
            self.exit()
            return
        
        if cmd in ["postest"]:
            self.log.info("Testing Source Positioning")
            self.positionTests()
            self.exit()
            return
        
        # Log Message generated by subfunction which understands how many lenslets
        self.generateAllLenslets()
        
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
    
    def exit(self):
        """Functions to run while shutting down the runner"""
        self.log.debug("Total Simulation took %2.3gs for %d lenslets with caching %s" % (time.clock() - self.start,
                len(self.lenslets),"enabled" if self.options.cache else "disabled"))
        self.log.info("Runner is done")
    
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
        import SEDInstrument, SEDSource
        with open(self.config["System"]["Configs"]["This"].rstrip(".yaml")+".dump.yaml",'w') as stream:
            yaml.dump(self.config["System"],stream,default_flow_style=False)
        
        Model = SEDInstrument.Instrument(self.config)
        Model.setup()
        Model.dumpConfig()
        
        Source = SEDSource.Source(self.config)
        Source.setup()
        Source.dumpConfig()
        
        
    def setup(self):
        """Performs all setup options"""
        self.setupModel()
        self.setupLenslets()
    
    
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
        import SEDSource
        
        self.Source = SEDSource.Source(self.config)
        
        self.Source.setup()
        
        self.Spectrum = self.Source.GetSpectrum()
        
        self.log.debug("Set Spectrum to %s" % self.Spectrum)
        
        
    
    def setupModel(self):
        """Sets up the SED Module Model"""
        if self.debug:
            start = time.clock()
        
        import SEDInstrument
        self.Limits = SEDInstrument.SEDLimits
        self.Model = SEDInstrument.Instrument(self.config)
        
        self.Model.plot = self.options.plot
        self.Model.setup()
        
        if self.debug:
            import matplotlib.pyplot as plt
            self.plt = plt
            end = time.clock()
            dur = end - start
            msg = "Setup took %1.5gs with caches %s." % (dur,"enabled" if self.options.cache else "disabled")
            self.log.debug(msg)
    
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
        handle = file(self.config["System"]["Dirs"]["Partials"] + "Positions.dat",'w')
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
            handle = file(self.config["System"]["Dirs"]["Partials"] + "Positions.dat",'a')
            np.savetxt(handle,np.array([np.hstack((corner,ncorner))]),fmt='%6.1f')
        except self.Limits:
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
            self.Model.cache_sed_subimage(i,spectrum,write=self.config["System"]["Cache"])
            
            if self.debug:
                self.Model.show()
                self.plt.savefig(self.config["System"]["Dirs"]["Partials"] + "%04d-Subimage.pdf" % i)
                self.plt.clf()
        
        except self.Limits:
            msg = "Skipped Lenslet %4d, Limits Out of Bounds" % i
            self.Model.log.info(msg)
        else:
            msg = "Cached Lenslet %4d" % i
            self.Model.log.info(msg)
        finally:
            self.prog.value += 1.0
        if self.debug:
            self.log.info("=>Finished Generating Spectrum %d" % self.prog.value)
        else:
            self.bar.render(int((self.prog.value/self.total) * 100),"L:%4d" % i)
        
    
    
    def generateAllLenslets(self):
        """Generate all lenslet spectra"""
        self.log.info("Generating Spectra in %d Lenslets" % len(self.lenslets))
        if not self.debug:
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
            self.Model.place_cached_sed(i,"Included Spectrum %d" % i,"Final",fromfile=self.config["System"]["Cache"])
            if self.debug:
                self.Model.show()
                self.plt.savefig(self.config["System"]["Dirs"]["Partials"] + "%04d-Fullimage.pdf" % i)
                self.plt.clf()
        except self.Limits:
            msg = "Encoutered Spectrum outside image boundaries %d" % i
            self.Model.log.info(msg)
        else:
            msg = "Placed Spectrum %d into image" % i
            self.Model.log.info(msg)
        finally:
            self.prog.value += 1
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
    


if __name__ == '__main__':
    Sim = Simulator()
    Sim.run()


