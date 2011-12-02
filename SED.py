#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SED.py
#  Simulation Software
#
#  Created by Alexander Rudy on 2011-11-04.
#  Copyright 2011 Alexander Rudy. All rights reserved.
#  Version 0.1.1
#


import numpy as np
import pyfits as pf
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.signal
import scipy.interpolate
import yaml

import os,logging,time,copy,collections

import logging.handlers

import AstroObject
from AstroObject.AstroSpectra import SpectraObject
from AstroObject.AstroImage import ImageObject,ImageFrame
from AstroObject.AnalyticSpectra import BlackBodySpectrum, AnalyticSpectrum, FlatSpectrum
from AstroObject.Utilities import *

from Utilities import *

__version__ = open(os.path.abspath(os.path.join(os.path.dirname(__file__),"VERSION")),'r').read()
__all__ = ["SEDLimits","Instrument"]


class SEDLimits(Exception):
    """A Basic Error-Differentiation Class.
    This error is used to express the fact that the SEDModel has encountered a spectrum which can't be placed as some part of it falls outside of the limits of the SED system."""
    pass


class Instrument(ImageObject):
    """This is a model container for the SEDMachine data simulator. This class is based on `AstroObject.ImageObject`, which it uses to provide some awareness of the way images are stored and retrieved. As such, this object has .save(), .select() and .show() etc. methods which can be used to examine the underlying data. It also means that this ImageObject subclass will contain the final, simulated image when the system is done.
    
    .. Note:: This object is not yet thread-safe, but coming versions of AstroObject will allow objects like this one to be locking and thread safe. Unfortunately, this limitation means that the simulator can generally not be run in multi-thread mode yet.
    
    """
    
    def __init__(self,config):
        super(Instrument, self).__init__()
        self.cache = False
        self.plot = True
        self.debug = False
        self.config = config
        self.defaults = []
        self.regenerate = not self.cache
        self.WLS = {}
        self.initLog()
        self.configure()
        self.setupLog()
    
    def update(self, d, u):
        """A deep update command for dictionaries.
        This is because the normal dictionary.update() command does not handle nested dictionaries."""
        for k, v in u.iteritems():
            if isinstance(v, collections.Mapping):
                r = self.update(d.get(k, {}), v)
                d[k] = r
            else:
                d[k] = u[k]
        return d
    
    def initLog(self):
        """Initializes the system logger. This logger starts with only a buffer, no actual logging output. The buffer is used to hold log messages before a logging output location has been specified."""
        self.log = logging.getLogger(__name__)
        
        if self.debug:
            self.logLevel = logging.DEBUG
        else:
            self.logLevel = logging.INFO
            
        self.log.setLevel(self.logLevel)
        logging.captureWarnings(True)
        self.filebuffer = logging.handlers.MemoryHandler(1e6) #Our handler will only handle 1-million messages... lets not go that far
        self.log.addHandler(self.filebuffer)
        self.consolebuffer = logging.handlers.MemoryHandler(1e6) #Our handler will only handle 1-million messages... lets not go that far
        self.log.addHandler(self.consolebuffer)
        
        self.log.info("--------------------------------")
        self.log.info("Welcome to the SED Machine model")
        self.log.debug(" Version %s" % __version__ )
    
    def setupLog(self):
        """Setup Logging Functions for the SEDMachine Model.
        
        This configures the logging system, including a possible console and file log. It then reads the logging buffer into the logfile.
        """
        
        if self.debug:
            self.logLevel = logging.DEBUG
        self.log.setLevel(self.logLevel)
        
        # Setup the Console Log Handler
        self.console = logging.StreamHandler()
        consoleFormat = self.config["Instrument"]["logging"]["console"]["format"]
        if self.config["Instrument"]["logging"]["console"]["level"]:
            self.console.setLevel(self.config["Instrument"]["logging"]["console"]["level"])
        else:
            self.console.setLevel(logging.INFO)
        consoleFormatter = logging.Formatter(consoleFormat)
        self.console.setFormatter(consoleFormatter)
        if self.config["Instrument"]["logging"]["console"]["enable"]:
            self.log.addHandler(self.console)
            self.consolebuffer.setTarget(self.console)
            self.consolebuffer.close()
            self.log.removeHandler(self.consolebuffer)
        
        
        
        dateFormat = "%Y-%m-%d-%H:%M:%S"
        self.logfile = None
        # Only set up the file log handler if we can actually access the folder
        if os.access(self.config["System"]["Dirs"]["Logs"],os.F_OK) and self.config["Instrument"]["logging"]["file"]["enable"]:
            filename = self.config["System"]["Dirs"]["Logs"] + self.config["Instrument"]["logging"]["file"]["filename"]+".log"
            self.logfile = logging.FileHandler(filename=filename,mode="a")
            self.logfile.setLevel(logging.DEBUG)
            fileformatter = logging.Formatter(self.config["Instrument"]["logging"]["file"]["format"],datefmt=dateFormat)
            self.logfile.setFormatter(fileformatter)
            self.log.addHandler(self.logfile)
            # Finally, we should flush the old buffers
            self.filebuffer.setTarget(self.logfile)
            self.filebuffer.close()
            self.log.removeHandler(self.filebuffer)
        
        
        self.log.debug("Configured Logging")
    
    
    # CONFIGURATION
    def configure(self):
        """This is a wrapper function for a variety of private configuration steps. The configuration generall happnes in stages:
        
        - Default Values are initalized for the instrument configuration
        - System loads the configuration file, updating the default values with the user-selected values. If no configuration file exists, the system will move on.
        - System loads the configuration from the Caches directory if caching is enabled.
        - System adds dynamic variables to the configuration. This is the feature which handles variable units etc.
        
        To see what the default configuration file looks like, run the runner script with the --dump argument. The runner script will requrie a subcommand, but dump will cause the program to exit before recieving any subcommands for action.
        
        You can force the script to ignore cached files in the runner script using the option `--no-cache`. To regenrate the cache manually, simply delete the contents of the Caches directory"""
        self._configureDefaults()
        self._configureFile()
        self._configureCaches()
        self._configureDynamic()
    
    def _configureDefaults(self):
        """Set up the default configure variable. If you change the default configuration variables in this function (instead of using a configuration file), the script will generally not detect the change, and so will not regenerate Cached files. You can force the script to ignore cached files in the runner script using the option `--no-cache`. To regenrate the cache manually, simply delete the contents of the Caches directory"""
        
        # Configuration Variables for The System
        self.config["Instrument"] = {}
        # Unit Conversion Defaults
        self.config["Instrument"]["convert"] = {}
        self.config["Instrument"]["convert"]["pxtomm"] = 0.0135
        # CCD / Image Plane Information
        self.config["Instrument"]["ccd_size"] = {}
        self.config["Instrument"]["ccd_size"]["px"] = 2048 #pixels
        self.config["Instrument"]["image_size"] = {}
        self.config["Instrument"]["image_size"]["mm"] = 40.0 #mm
        # Telescope Information
        self.config["Instrument"]["tel_radii"] = {}
        self.config["Instrument"]["tel_radii"]["px"] = 2.4 / 2.0
        self.config["Instrument"]["tel_obsc"] = {}
        self.config["Instrument"]["tel_obsc"]["px"] = 0.4 / 2.0
        # PSF Information
        self.config["Instrument"]["psf_size"] = {}
        self.config["Instrument"]["psf_size"]["px"] = 0
        # For a gaussian PSF
        self.config["Instrument"]["psf_stdev"] = {}
        self.config["Instrument"]["psf_stdev"]["px"] = 1.0
        # Image Generation Density
        self.config["Instrument"]["density"] = 5
        self.config["Instrument"]["padding"] = 5
        # Default Gain Value
        self.config["Instrument"]["gain"] = 1e-6
        # Noise Information
        self.config["Instrument"]["dark"] = 20 # counts per pixel per second at some fixed degree c
        self.config["Instrument"]["bias"] = 20 # counts per pixel at some fixed degree c
        self.config["Instrument"]["exposure"] = 120 #Seconds
        # File Information for data and Caches
        self.config["Instrument"]["files"] = {}
        self.config["Instrument"]["files"]["lenslets"] = "Data/xy_17nov2011_v57.TXT"
        self.config["Instrument"]["files"]["dispersion"] = "Data/dispersion_12-10-2011.txt"
        self.config["Instrument"]["files"]["encircledenergy"] = "Data/encircled_energy_4nov11.TXT"
        # MPL Plotting Save Format
        self.config["Instrument"]["plot_format"] = ".pdf"
        
        # Logging Configuration
        self.config["Instrument"]["logging"] = {}
        self.config["Instrument"]["logging"]["console"] = {}
        self.config["Instrument"]["logging"]["console"]["enable"] = False
        self.config["Instrument"]["logging"]["console"]["format"] = "... ...%(message)s"
        self.config["Instrument"]["logging"]["console"]["level"] = False
        self.config["Instrument"]["logging"]["file"] = {}
        self.config["Instrument"]["logging"]["file"]["enable"] = True
        self.config["Instrument"]["logging"]["file"]["filename"] = "SED"+"-"+time.strftime("%Y-%m-%d")
        self.config["Instrument"]["logging"]["file"]["format"] = "%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s"
        
        self.defaults += [copy.deepcopy(self.config)]
        
        self.log.debug("Set Instrument Configuration Defaults")
    
    def _configureFile(self):
        """Attempt to load the default configuration from the working directory"""
        FileName = self.config["System"]["Configs"]["Instrument"]
        try:
            stream = open(FileName,'r')
        except IOError:
            self.log.warning("Configuration File Not Found: %s" % FileName)
        else:
            self.update(self.config["Instrument"],yaml.load(stream))
            stream.close()
            self.log.debug("Loaded Configuration from %s" % FileName)
        finally:
            self.defaults += [copy.deepcopy(self.config)]
    
    def _configureCaches(self):
        """If we are using the caching system, set up the configuration cache.
        
        First, we set default Caching variables. Next, load these same variables (and overwrite) from the system configuration. Next, we load the cached configuration. If the cached configuration updates the generated configuration in any way, """
        if not self.cache:
            self.log.debug("Skipping Cache Configuration")
            return
        
        # Default Cache Variables
        CacheFiles = {}
        CacheFileBase = "SED.Instrument"
        CacheFiles["telescope"] = self.config["System"]["Dirs"]["Caches"] + CacheFileBase + ".tel"    + ".npy"
        CacheFiles["psf"]       = self.config["System"]["Dirs"]["Caches"] + CacheFileBase + ".psf"    + ".npy"
        CacheFiles["conv"]      = self.config["System"]["Dirs"]["Caches"] + CacheFileBase + ".conv"   + ".npy"
        CacheFiles["config"]    = self.config["System"]["Dirs"]["Caches"] + CacheFileBase + ".config" + ".yaml"
        CacheFiles["wls"]       = self.config["System"]["Dirs"]["Caches"] + CacheFileBase + ".wls"    + ".npy"
        
        CacheFiles = self.update(CacheFiles,self.config["System"]["CacheFiles"])
        self.config["System"]["CacheFiles"] = CacheFiles
        
        self.log.debug("Configured Cache Variables")
        self.defaults += [copy.deepcopy(self.config)]
        
        # Load the Cached Configuration
        FileName = self.config["System"]["CacheFiles"]["config"]
        try:
            stream = file(FileName,'r')
            self.update(self.config["Instrument"],yaml.load(stream))
        except IOError:
            self.log.info("Cached Configuration File Not Found: %s" % FileName)
        except TypeError:
            self.log.critical("Cached Configuration File has a Problem... skipping.")
            stream.close()
        else:
            stream.close()
        
        # If the cache is different from the generated configuration, then don't trust any caches
        if self.defaults[-1]["Instrument"] == self.config["Instrument"]:
            self.regenearte = False
            self.log.debug("Configuration has not changed, will not regenerate data.")
        else:
            self.regenearte = True
            self.config["Instrument"] = self.defaults[-1]["Instrument"]
            self.log.info("Configuration appears to have changed, will regenerate data.")
        
        return
    
    def _configureDynamic(self):
        """Generate dynamic configuration values.
        
        Currently, the configuration variables that use \"px\" or \"mm\" keys automatically have their counter-part filled in. As well, the conversion has its counterpart filled in. To disable this dynamic value setting, include a \"calc:False\" variable in the configuration.
        
        Example Configuration::
        
            image_size:
                mm: 40 # system will calculate pixels
            ccd_size:
                px: 2048
                calc: False # Conversion will not be performed by system
                mm: 30
        
        
        """
        self.defaults += [copy.deepcopy(self.config)]
        
        if "calc" not in self.config["Instrument"]["convert"]:
            if "mmtopx" not in self.config["Instrument"]["convert"] and "pxtomm" in self.config["Instrument"]["convert"]:
                self.config["Instrument"]["convert"]["mmtopx"] = 1.0 / self.config["Instrument"]["convert"]["pxtomm"]
                self.config["Instrument"]["convert"]["calc"] = True
            else:
                self.config["Instrument"]["convert"]["pxtomm"] = 1.0 / self.config["Instrument"]["convert"]["mmtopx"]
                self.config["Instrument"]["convert"]["calc"] = True
        
        self.config["Instrument"] = self._setUnits(self.config["Instrument"],None)
        
        self.config["Instrument"]["image_size"]["px"] = np.round( self.config["Instrument"]["image_size"]["px"] , 0 )
    
    def _setUnits(self,config,parent):
        """docstring for _setUnits"""
        r = copy.deepcopy(config)
        for k, v in config.iteritems():
            if isinstance(v, collections.Mapping):
                r[k] = self._setUnits(r[k],k)
            elif k == "mm":
                if ("calc" in r) and ("px" in r):
                    pass
                elif ("calc" not in config):
                    r["px"] = v * self.config["Instrument"]["convert"]["mmtopx"]
                    r["calc"] = True
                else:
                    self.log.warning("Value for %s set in both px and mm." % parent)
            elif k == "px":
                if ("calc" in r) and ("mm" in r):
                    pass
                elif ("calc" not in r):
                    r["mm"] = v * self.config["Instrument"]["convert"]["pxtomm"]
                    r["calc"] = True
                else:
                    self.log.warning("Value for %s set in both px and mm." % parent)
        return r
    
    
    # Cacheing Functions
    def regenerateCache(self):
        """Cache calculated components of the system, including the telescope image and encircled energy image. Caches are stored to speed up system initalization. This function regenerates all cached files, including the configuration file. You can force the script to ignore cached files in the runner script using the option `--no-cache`. To regenrate the cache manually, simply delete the contents of the Caches directory"""
        with open(self.config["System"]["CacheFiles"]["config"],'w') as stream:
            yaml.dump(self.defaults[-1]["Instrument"],stream,default_flow_style=False)
        
        np.save(self.config["System"]["CacheFiles"]["telescope"],self.TELIMG)
        np.save(self.config["System"]["CacheFiles"]["psf"],self.PSFIMG)
        np.save(self.config["System"]["CacheFiles"]["conv"],self.FINIMG)
    
    def cachedKernel(self):
        """Load cached kernels from the Caches directory. If any file is missing, it will attempt to trigger regeneration of the cache.
        You can force the script to ignore cached files in the runner script using the option `--no-cache`. To regenrate the cache manually, simply delete the contents of the Caches directory"""
        try:
            # Telescope Image Setup
            self.TELIMG = np.load(self.config["System"]["CacheFiles"]["telescope"])
            # PSF Setup
            self.PSFIMG = np.load(self.config["System"]["CacheFiles"]["psf"])
            # Preconvolved System
            self.FINIMG = np.load(self.config["System"]["CacheFiles"]["conv"])
        except IOError as e:
            self.log.warning("Cached files not found, using configuration to generate files. Error: %s" % e )
            self.regenerate = True
        else:
            self.log.debug("Loaded Telescope Images for Numpy Files")
        
    
    def cachedWL(self):
        """Load cached wavelengths from the Caches directory. If any file is missing, it will attempt to trigger regeneration of the cache.
        You can force the script to ignore cached files in the runner script using the option `--no-cache`. To regenrate the cache manually, simply delete the contents of the Caches directory"""
        try:
            # Cached Wavelengths
            cachedWL = np.load(self.config["System"]["CacheFiles"]["wls"]+".npy")
            self.WLS = dict(cachedWL)
        except IOError as e:
            self.log.warning("Cached files not found, using configuration to generate files. Error: %s" % e )
            self.regenerate = True
        else:
            self.log.debug("Loaded Wavelengths from Numpy Files")
    
    def resetWLCache(self):
        """Resets the Wavelength Cache, which will force it to regenerate during the simulation"""
        self.log.debug("Forcing Wavelenghts to Regenerate")
        self.WLS = {}
    
    def dumpConfig(self):
        """Dumps a valid configuration file on top of any old configuration files. This is useful for examining the default configuration fo the system, and providing modifications through this interface."""
        with open(self.config["System"]["Configs"]["Instrument"].rstrip(".yaml")+".dump.yaml",'w') as stream:
            yaml.dump(self.defaults[1]["Instrument"],stream,default_flow_style=False)
    
    def setup(self):
        """After the object has been initialized, it must be setup. Setup relies on an established configuration to determine what fixed parts of the system should be generated. Actions taken in the setup phase are:
        
        - Load data from the Cache (if enabled)
        - Regenerate Kernel Functions (if required)
        - Regenerate Cache Files (if required, and if enabled)
        - Load optical system data (not cached right now)
        - Generate a blank image (no need to cache... np.zeros is fast!)
        
        """
        
        if self.cache:
            self.cachedKernel()
            self.cachedWL()
        if self.regenerate:
            self.regenerateKernel()
            self.resetWLCache()
        if self.cache and self.regenerate:
            self.regenerateCache()
        
        if self.log.getEffectiveLevel() <= logging.DEBUG and self.plot:
            self.plotKernelPartials()
        
        # Get layout data from files
        self.loadOpticsData(self.config["Instrument"]["files"]["lenslets"],self.config["Instrument"]["files"]["dispersion"])
        # TODO: We could cache this...
        
        self.generate_blank()
        
        self.log.info("Done with SEDM setup")
    
    
    def plotKernalPartials(self):
        """Plots the kernel data partials"""
        self.log.debug("Generating Kernel Images")
        plt.clf()
        plt.imshow(self.TELIMG)
        plt.title("Telescope Image")
        plt.colorbar()
        plt.savefig("%sTelImage%s" % (self.config["System"]["Dirs"]["Partials"],self.config["Instrument"]["plot_format"]))
        plt.clf()
        plt.imshow(self.PSFIMG)
        plt.title("PSF Image")
        plt.colorbar()
        plt.savefig("%sPSFImage%s" % (self.config["System"]["Dirs"]["Partials"],self.config["Instrument"]["plot_format"]))
        plt.clf()
        plt.imshow(self.FINIMG)
        plt.title("Convolved Tel + PSF Image")
        plt.colorbar()
        plt.savefig("%sFINImage%s" % (self.config["System"]["Dirs"]["Partials"],self.config["Instrument"]["plot_format"]))
        plt.clf()
    
    # Kernel Creation for Image Manipulation
    def regenerateKernel(self):
        """Regenerate a kernel from the configuration values."""
        self.log.info("Generating Kernels for Image System")
        # Telescope Image Setup
        self.TELIMG = self.get_tel_kern()
        # PSF Setup
        self.PSFIMG = self.get_psf_kern()
        # Preconvolved System
        self.FINIMG = sp.signal.convolve(self.PSFIMG,self.TELIMG,mode='same')
    
    def psf_kern(self,filename,size=0,truncate=False,header_lines=18):
        """Generates a PSF Kernel from a file with micron-encircled energy conversions. The file should have two columns, first, microns from the center of the PSF, and second, the fraction of encircled energy at that distance from the PSF.
        
        The calculation is then carried out using a spline fit to this data. From the spline fit, the function returns the first derivative of the encircled energy at each point. This in effect is the amount of energy at each point. These values are then normalized, to create a PSF mask for the instrument.
        
        The `size` parameter specifies the size of the kernel to use. If the size is greater than the encircled energy data, then a larger figure will be returned. If `size` is smaller than the encircled energy data, it will return an image the size of the encircled energy data, unless the `truncate` parameter is set to `true`.
        
        The `header_lines` parameter defaults to 18, which works with Zemax encircled energy output."""
        
        uM,FR = np.genfromtxt(filename,skip_header=header_lines).T
        # Convert microns to milimeters, then pixels, then dense pixels
        PX = uM * 1e-3 * self.config["Instrument"]["convert"]["mmtopx"] * self.config["Instrument"]["density"]
        # Set the frame size for the PSF
        if np.max(PX) <= size or truncate:
            size = np.int(size)
        else:
            size = np.int(np.max(PX))
        # Create the Interpolation Function
        fit_vars = sp.interpolate.splrep(PX,FR)
        fit = lambda x : sp.interpolate.splev(x,fit_vars,der=1)
        vfit = np.vectorize(fit)
        # Create the 2-D function application grid
        x , y = np.mgrid[-size:size+1,-size:size+1]
        # Convert this grid into the distance from the center of the PSF at each point
        r = np.sqrt(x**2 + y**2)
        v = vfit(r)
        val = v
        self.log.debug("Generated a PSF Kernel for the encircled energy file %s with shape %s" % (filename,str(v.shape)))
        return val / np.sum(val)
    
    def circle_kern(self,radius,size=0,normalize=False):
        """Generate a Circle Kernel for modeling the \"Image of the Telescope\". The radius should be set in array units.
        
        `size` will determine the size of the array image, unless `size` is less than `radius`, in which case the image will be automatically increased to fit the entire circle.
        
        `normalize` controls whether the data is normalized or not. If it is not normalized, the data will have only 1.0 and 0.0 values, where 1.0 is within the radius, and 0.0 is outside the raidus."""
        if size < radius:
            size = int(radius)
        else:
            size = int(size)
        radius = int(radius)
        x, y = np.mgrid[-size:size+1, -size:size+1]
        d = np.sqrt(x**2.0 + y**2.0)
        v = (d <= radius).astype(np.float)
        if normalize:
            return v / np.sum(v)
        else:
            return v
    
    def gauss_kern(self,stdev,size=0,stdevy=None,sizey=0):
        """ Returns a normalized 2D gaussian kernel array for convolutions.
        
        `stdev` is the standard deviation in the x-direction. If the `stdevy` keyword is not set, then it will be used as the standard deviation in the y-direction as well.
        
        `size` will determine the size of the returned image unless `size` is less than `stdev**2`.
        
        Results from this function are always normalized.
        """
        if size < (stdev**2.0):
            size = np.int(stdev**2.0)
        else:
            size = np.int(size)
        if not stdevy:
            stdevy = stdev
        if sizey < (stdevy**2.0):
            sizey = np.int(stdevy**2.0)
        else:
            sizey = np.int(sizey)
        
        x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
        g = np.exp(-(x**2.0/np.float(stdev**2.0)+y**2.0/np.float(stdevy**2.0)))
        
        return g / g.sum()
    
    
    
    def get_tel_kern(self):
        """Returns the telescope kernel. This kernel is built by creating a circle mask for the size of the telescope mirror, and then subtracting a telescope obscuration from the center of the mirror image. The values for all of these items are set in the configuration file."""
        TELIMG = self.circle_kern( self.config["Instrument"]["tel_radii"]["px"] * self.config["Instrument"]["density"] )
        center = self.circle_kern( self.config["Instrument"]["tel_obsc"]["px"] * self.config["Instrument"]["density"] ,
            self.config["Instrument"]["tel_radii"]["px"] * self.config["Instrument"]["density"] , False )
        TELIMG -= center
        TELIMG = TELIMG / np.sum(TELIMG)
        self.log.debug("Generated a Telescpe Kernel with shape %s" % (str(TELIMG.shape)))
        return TELIMG
    
    def get_psf_kern(self):
        """Returns the PSF Kernel. The function first tries to read the encircled energy file. In this case, if a `psf_size` is set in the instrument configuration, this value will be used to truncate the size of the encircled energy function. If the encircled energy function cannot be loaded, the system will fall back on to a gaussian psf as configured by the instrument."""
        if self.config["Instrument"]["psf_size"]["px"] != 0:
            size = self.config["Instrument"]["psf_size"]["px"] * self.config["Instrument"]["density"]
            truncate = True
        else:
            size = 0
            truncate = False
        try:
            PSFIMG = self.psf_kern( self.config["Instrument"]["files"]["encircledenergy"],size,truncate)
        except IOError as e:
            self.log.warning("Could not access encircled energy file: %s" % e)
            PSFIMG = self.gauss_kern( (self.config["Instrument"]["psf_stdev"]["px"] * self.config["Instrument"]["density"]) )
        else:
            self.log.debug("Loaded Encircled Energy from %s" % self.config["Instrument"]["files"]["encircledenergy"])
        return PSFIMG
    
    
    def get_blank_img(self):
        """Returns an image of the correct size to use for spectrum placement. Set by the `image_size` value."""
        return np.zeros((self.config["Instrument"]["image_size"]["px"],self.config["Instrument"]["image_size"]["px"]))
    
    
    # Wavelength Functions
    def get_wavelengths(self,lenslet_num):
        """Returns a tuple of ([x,y],wavelength,delta wavelength). Thus, this function calculates the x,y pixel positions of the spectra, the wavelength at each pixel, and the delta wavelength for each pixel.
        
        .. Note:: This function can use cached calculations, which do not seem to necessarily speed up calculation.
        
        .. Note:: The logic of this collection is in `_get_wavelengths()`, which is heavily documented in source code.
        """
        if lenslet_num in self.WLS:
            self.log.debug("Using Cached WLs for %d" % lenslet_num)
            results = self.WLS[lenslet_num]
        else:
            self.log.debug("Regenerating WLs for %d" % lenslet_num)
            results = self._get_wavelengths(lenslet_num)
            self.WLS[lenslet_num] = np.array(results)
            self.regenerate = True
        return results
    
    def positionCaching(self):
        """Save the wavelength and position cache. This must be done after all wavelengths have been collected."""
        if self.cache and self.regenerate:
            toCache = np.array(zip(self.WLS.keys(),self.WLS.values()))
            np.save(self.config["System"]["CacheFiles"]["wls"],toCache)
        
    
    def _get_wavelengths(self,lenslet_num):
        """Returns an array of wavelengths for a given lenslet number. See `get_wavelenghts()`.
        
        .. Note:: This function is heavily documented in the source code.
        """
        # This method will be the critical slow point for the whole system...
        # It should be re-written to do the conversions in a better way, but
        # I'm not sure how to do that right now.
        
        # First, we take only data points which apply to this lenslet
        use = lenslet_num == self.ix
        
        # Each lenslet should have three points. It might not though, because we clipped some lenselet points that were too close to the edge.
        # We really should move this logic up higher.
        if len(self.xpix[use]) < 3:
            raise SEDLimits
        
        if np.any(self.xpix[use] == 0):
            raise SEDLimits
        
        if np.any(self.ypix[use] == 0):
            raise SEDLimits
        
        # We ignore anything that has an x-dispersion across more than 30 pixels.
        if np.any(np.abs(np.diff(self.xpix[use])) > 30):
            raise SEDLimits
        
        # Interpolation to convert from wavelength to pixels.
        #   The accuracy of this interpolation is not important.
        #   Rather, it is used to find the pixels where the light will fall
        #   and is fed an array that is very dense, used on this dense interpolation
        #   and then binned back onto pixels. Thus it will be used to get a list
        #   of all illuminated pixels.
        fx = np.poly1d(np.polyfit(self.lams[use], self.xpix[use], 2))
        fy = np.poly1d(np.polyfit(self.lams[use], self.ypix[use], 2))
        
        # Find the starting and ending position of the spectra
        startix = np.argmin(self.lams[use])
        endix = np.argmax(self.lams[use])
        start = np.array([self.xs[use][startix],self.ys[use][startix]])
        end = np.array([self.xs[use][endix],self.ys[use][endix]])
        
        # Get the total length of the spectra
        distance = np.sqrt(np.sum(end-start)**2)
        
        if distance == 0:
            raise SEDLimits
        
        # Find the length in units of (int) pixels
        npix = (distance * self.config["Instrument"]["convert"]["mmtopx"]).astype(np.int) * self.config["Instrument"]["density"]
        
        # Create a data array one hundred times as dense as the number of pixels
        #   This is the super dense array which will use the above interpolation
        superDense_lam = np.linspace(np.min(self.lams[use]),np.max(self.lams[use]),npix*100)
        
        # Interpolate along our really dense set of wavelengths to find all possible
        # illuminated pixel positions in this spectrum
        superDense_pts = np.array([fx(superDense_lam),fy(superDense_lam)])
        
        # Measure the distance along our really dense set of points
        superDense_interval = np.sqrt(np.sum(np.power(np.diff(superDense_pts,axis=1),2),axis=0))
        superDense_distance = np.cumsum(superDense_interval)
        
        # Adjust the density of our points. This rounds all values to only full pixel values.
        if self.config["Instrument"]["density"] == 1:
            superDense_pts = superDense_pts.astype(np.int)
        else:
            superDense_pts = np.round(superDense_pts * self.config["Instrument"]["density"]) / self.config["Instrument"]["density"]
            superDense_int = (superDense_pts * self.config["Instrument"]["density"]).astype(np.int)
        
        # We can identify unique points using the points when the integer position ratchets up or down.
        unique_x,unique_y = np.diff(superDense_int).astype(np.bool)
        
        # We want unique index to include points where either 'y' or 'x' ratchets up or down
        unique_idx = np.logical_or(unique_x,unique_y)
        
        # Remove any duplicate points. This does not do so in order, so we must
        # sort the array of unique points afterwards...
        unique_pts = superDense_pts[:,1:][:,unique_idx]
        
        # An array of distances to the origin of this spectrum, can be used to find wavelength
        # of light at each distance
        distance = superDense_distance[unique_idx] * self.config["Instrument"]["convert"]["pxtomm"]
        
        # Re sort everything by distnace along the trace.
        # Strictly, this shouldn't be necessary if all of the above functions preserved order.
        sorted_idx = np.argsort(distance)
        
        # Pull out sorted valuses
        distance = distance[sorted_idx]
        points = unique_pts[:,sorted_idx].T
        self.log.debug("Points set using original, superDense array, bounds [%1.4f,%1.4f]"
            % (np.min(points),np.max(points)))
        
        
        # Pull out the original wavelengths
        wl_orig = superDense_lam[unique_idx][sorted_idx]
        wl = wl_orig
        self.log.debug("Wavelengths set using original, superDense array, bounds [%1.4f,%1.4f]"
            % (np.min(wl),np.max(wl)))
        # We are getting some odd behavior, where the dispersion function seems to not cover the whole
        # arc length and instead covers only part of it. This causes much of our arc to leave the desired
        # and available wavelength boundaries. As such, I'm disabling the more accurate dispersion mode.
        
        # Convert to wavelength space along the dispersion spline.
        # wl = self.spline(distance)
        
        # This is a debugging area which will make a lot of graphs and annoying messages.
        if self.log.getEffectiveLevel() <= logging.DEBUG and self.plot:
            self.log.debug("Plotting Wavelength Information")
            # This graph shows the change in distance along arc per pixel.
            # The graph should produce all data points close to each other, except a variety of much lower
            # data points which are caused by the arc crossing between pixels.
            plt.clf()
            plt.title("$\Delta$Distance Along Arc")
            plt.xlabel("x (px)")
            plt.ylabel("$\Delta$Distance along arc (px)")
            plt.plot(points[:-1,1],np.diff(distance) * self.config["Instrument"]["convert"]["mmtopx"],'g.')
            plt.savefig("%s%04d-Distances_Diff%s" % (self.config["System"]["Dirs"]["Partials"],lenslet_num,self.config["Instrument"]["plot_format"]))
            plt.clf()
        
        return points,wl,np.diff(wl)
    
    
    # Data Loading Functions
    def loadDispersionData(self,dispspec):
        """This loads the dispersion data. Dispersion data is provided in mm along the spectrum and wavelength. This method saves both the raw data (in :attr:`self.MM` and :attr:`self.WL`) as well as functions which can be used to convert across the dispersion. Functions by default take a mm position and return the wavelength at this position."""
        WLtoMM = np.genfromtxt(dispspec).T
        WL = WLtoMM[0]
        MM = WLtoMM[1]
        MM -= np.min(MM) #Zero offset the MM data
        self.MM,self.WL = MM,WL
        
        # A note about functions:
        #
        # All of these functions will return the wavelength of a position, in milimeters,
        # from the origin of the spectrum. We define the origin of the spectrum as the
        # postion where the 0.37 micron wavelength falls.
        
        # Define a basic Floor function interpolation
        self.floor = lambda x: np.vectorize(lambda ax: MM[(np.abs(ax-WL)).argmin()])(x)
        # Define a piecewise linear interpolation
        self.linear = sp.interpolate.interp1d(WL,MM)
        # Generate the SPLINE function:
        self.spvars = sp.interpolate.splrep(MM,WL)
        self.speval = sp.interpolate.splev
        self.speval.bounds_error = True
        self.spline = lambda x: self.speval(x,self.spvars)
    
    def loadLensletData(self,laspec):
        """This function loads data about lenslet positions, and thier dispersion through the prism. The data are original produced by Zeemax. This function reads the Zeemax data directly and then cleans the data in certain ways, preparing it for use later in the system.
        
        ..Note:: The source of this function is well documented.
        
        ..Note:: This function does not store variables neatly. As such, it has no built-in caching system.
        """
        # Load Lenslet Specification File
        ix, p1, p2, lams, xs, ys = np.genfromtxt(laspec,skip_header=1).T
        # This data describes the following:
        # ix - Index (number)
        # p1 - Pupil position in the x-direction
        # p2 - Pupil position in the y-direction
        # lams - wavelengths for this position
        # xs - X position (in mm, offest from top right corner) of this wavelength
        # ys - Y Position (in mm, offset from top right corner) of this wavelength
        
        # Correctly Type Lenslet Specification Data
        ix = ix.astype(np.int) #Indicies should always be integers
        
        # Center xs and ys on detector with a corner at 0,0
        xs += (self.config["Instrument"]["image_size"]["mm"]/2)
        ys += (self.config["Instrument"]["image_size"]["mm"]/2)
        
        # Find the xs and ys that are not within 0.1 mm of the edge of the detector...
        ok = (xs > 0.1) & (xs < self.config["Instrument"]["image_size"]["mm"]-0.1) & (ys > 0.1) & (ys < self.config["Instrument"]["image_size"]["mm"]-0.1)
        ix, p1, p2, lams, xs, ys = ix[ok], p1[ok], p2[ok], lams[ok], xs[ok], ys[ok]
        # We remove these positions because we don't want to generate spectra for them.
        
        # This simply generates a list of all of the lenslets
        self.lenslets = np.unique(ix)
        
        # Convert the xs and ys to pixel positions
        xpix = np.round(xs * self.config["Instrument"]["convert"]["mmtopx"],0).astype(np.int)
        ypix = np.round(ys * self.config["Instrument"]["convert"]["mmtopx"],0).astype(np.int)
        
        # Determine the center of the whole system by finding the x position that is closest to 0,0 in pupil position
        cntix = np.argmin(p1**2 + p2**2)
        self.center = (xs[cntix] * self.config["Instrument"]["convert"]["mmtopx"], ys[cntix] * self.config["Instrument"]["convert"]["mmtopx"])
        
        self.ix, self.p1, self.p2, self.lams, self.xs, self.ys = ix, p1, p2, lams, xs, ys
        self.xpix, self.ypix = xpix, ypix
    
    def loadOpticsData(self,laspec,dispspec):
        """Loads an optical conversion based on the lenslet array spec and dispersion spec files provided. This wrapper function handles both file loading functions. See `loadDispersionData` and `loadLensletData`."""
        self.loadDispersionData(dispspec)
        self.loadLensletData(laspec)
    
    
    # Image Tracing Function
    def get_dense_image(self,lenslet,spectrum):
        """This function returns a dense image array, with flux placed into single pixels. The image returned is too dense for use on the final CCD, and must be also be convolved with the telescope image and spectral PSF. The image, and the numpy array first corner is returned by this function.
        
        .. Note:: The source code of this function is heavily documented."""
        points, wl, deltawl = self.get_wavelengths(lenslet)
        # Some notes about this function:
        #  points:  These points are in normal CCD Pixel units, but are not integers,
        #           as they correspond to dense pixels on the overdense images
        #  wl:      A list of wavelengths that correspond to the points above
        #  deltawl: The np.diff(wl) result, which can be used to handle the flux per pixel
        
        # Create the spectra, by calling the spectra with wavelength
        # (in meters, note the conversion from wl, which was originally in microns)
        # Then use the deltawl to get the true amount of flux in each pixel
        self.log.debug("Asking for spectrum with bounds [%1.4e,%1.4e]" % (np.max(wl),np.min(wl)))
        radiance = spectrum(wl[:-1]*1e-6) * self.config["Instrument"]["gain"]
        self.log.debug("Generated spectrum with bounds [%1.4e,%1.4e]" % (np.max(radiance),np.min(radiance)))
        self.log.debug("Re-scaling Radiance by deltawl with bounds [%1.4e,%1.4e]" % (np.max(deltawl),np.min(deltawl)))
        flux = radiance[1,:] * deltawl
        self.log.debug("Got Flux with bounds [%1.4e,%1.4e]" % (np.max(flux),np.min(flux)))
        
        # This debugging area generates plots for us.
        if self.log.getEffectiveLevel() <= logging.DEBUG and self.plot:
            self.log.debug("Generating Plots for Spectra...")
            plt.clf()
            plt.plot(wl[:-1],flux,"b.")
            plt.title("Generated, Fluxed Spectra")
            plt.xlabel("Wavelength ($\mu m$)")
            plt.ylabel("Flux (Units undefined)")
            plt.savefig("%s%04d-SpecFlux%s" % (self.config["System"]["Dirs"]["Partials"],lenslet,self.config["Instrument"]["plot_format"]))
            plt.clf()
            plt.plot(wl[:-1],deltawl,"g.")
            plt.title("$\Delta\lambda$ for each pixel")
            plt.xlabel("Wavelength ($\mu m$)")
            plt.ylabel("$\Delta\lambda$ per pixel")
            plt.savefig("%s%04d-SpecDeltaWL%s" % (self.config["System"]["Dirs"]["Partials"],lenslet,self.config["Instrument"]["plot_format"]))
            plt.clf()
        
        # Take our points out. Note from the above that we multiply by the density in order to do this
        xorig,yorig = (points * self.config["Instrument"]["density"])[:-1].T.astype(np.int)
        x,y = (points * self.config["Instrument"]["density"])[:-1].T.astype(np.int)
        # Get the way in which those points correspond to actual pixels.
        # As such, this array of points should have duplicates
        xint,yint = points.T.astype(np.int)
        
        # Zero-adjust our x and y points. They will go into a fake subimage anyways, so we don't care
        # for now where they would be on the real image
        x -= np.min(x)
        y -= np.min(y)
        
        # Get the approximate size of our spectra
        xdist = np.max(x)-np.min(x)
        ydist = np.max(y)-np.min(y)
        
        # Convert this size into an integer number of pixels for our subimage. This makes it
        # *much* easier to register our sub-image to the master, larger pixel image
        xdist += (self.config["Instrument"]["density"] - xdist % self.config["Instrument"]["density"])
        ydist += (self.config["Instrument"]["density"] - ydist % self.config["Instrument"]["density"])
        
        # Move our x and y coordinates to the middle of our sub image by applying padding below each one.
        x += self.config["Instrument"]["padding"] * self.config["Instrument"]["density"]
        y += self.config["Instrument"]["padding"] * self.config["Instrument"]["density"]
        
        # Find the first (by the flatten method) corner of the subimage,
        # useful for placing the sub-image into the full image.
        corner = np.array([ xint[np.argmax(x)], yint[np.argmin(y)]])
        self.log.debug("Corner Position in Integer Space: %s" % corner)
        corner *= self.config["Instrument"]["density"]
        realcorner = np.array([ xorig[np.argmax(x)], yorig[np.argmin(y)]])
        offset = corner - realcorner
        corner /= self.config["Instrument"]["density"]
        self.log.debug("Corner Position Offset in Dense Space: %s" % (offset))
        if self.log.getEffectiveLevel() <= logging.DEBUG:
            with open("Partials/Offsets.dat",'a') as handle:
                np.savetxt(handle,offset.T)
        corner -= np.array([-self.config["Instrument"]["padding"],self.config["Instrument"]["padding"]])
        
        x += offset[0]
        y += offset[1]
        
        # Create our sub-image, using the x and y width of the spectrum, plus 2 padding widths.
        # Padding is specified in full-size pixels to ensure that the final image is an integer
        # number of full-size pixels across.
        xsize = xdist+2*self.config["Instrument"]["padding"]*self.config["Instrument"]["density"]
        ysize = ydist+2*self.config["Instrument"]["padding"]*self.config["Instrument"]["density"]
        img = np.zeros((xsize,ysize))
        
        # Place the spectrum into the sub-image
        img[x,y] = flux
        self.log.debug("Placing flux (shape: %s ) for spectrum %4d into a sub-image (shape: %s)."
            % (flux.shape,lenslet,img.shape))
        
        if self.log.getEffectiveLevel() <= logging.DEBUG:
            np.savetxt("Partials/SubimageValues.dat",np.array([x,y,wl[:-1],deltawl,flux]).T)
        
        return img, corner
    
    
    # Image Convolution
    def get_sub_image(self,lenslet,spectrum,fast=False):
        """Returns a sub-image for a given lenslet, which has been binned and convolved with the appropriate telescope and PSF images. There are two modes for this fucntion: fast, which uses a single convolution, and slow, which does each convolution separately and returns intermediate stage images for diagnostics."""
        # This function gets the dense sub-image with the spectrum placed in single dense pixels
        img,corner = self.get_dense_image(lenslet,spectrum)
        self.log.debug("Retrieved Dense Image for %4d" % lenslet)
        
        if fast:
            # Convolve with the PSF and Telescope Image simultaneously
            img2 = sp.signal.convolve(img,self.FINIMG,mode='same')
            self.log.debug("Convolved Dense Image with PSF and Telescope for %4d" % lenslet)
            # Bin the image back down to the final pixel size
            small = bin(img2,self.config["Instrument"]["density"]).astype(np.int16)
            self.log.debug("Binned Dense Image to Actual Size for %4d" % lenslet)
            return small,corner
        else:
            # Convolve this spectrum with an appropriate image of the telescope
            img_tel = sp.signal.convolve(img,self.TELIMG,mode='same')
            self.log.info("Convolved Dense Image with Telescope for %4d" % lenslet)
            # Convolve again with an appropriate PSF
            img_tel_psf = sp.signal.convolve(img_tel,self.PSFIMG,mode='same')
            self.log.info("Convolved Dense Image with PSF for %4d" % lenslet)
            # Bin the image back down to the final pixel size
            small = bin(img_tel_psf,self.config["Instrument"]["density"]).astype(np.int16)
            self.log.info("Binned Dense Image for %4d" % lenslet)
            return small,corner,(img,img_tel,img_tel_psf)
    
    
    # Image Caching
    def cache_sed_subimage(self,lenslet,spectrum,write=False,do_return=False):
        """Generates a sub image, and saves that result to this object. Should be thread-safe."""
        if self.log.getEffectiveLevel() <= logging.DEBUG:
            small, corner, steps = self.get_sub_image(lenslet,spectrum)
        else:
            small, corner = self.get_sub_image(lenslet,spectrum,fast=True)
        
        self.log.info("Retrieved SED Subimage for lenslet %d" % lenslet)
        
        label = "SUBIMG%d" % lenslet
        
        self.save(small,label)
        self.frame().header=dict(Lenslet=lenslet,Spectrum=spectrum.label)
        self.frame().metadata=dict(Corner=corner)
        
        Stages = ["Raw Image","Convolved with Telescope","Convolved with Telescope and PSF"]
        StagesF = ["Raw","Tel","PSF"]
        if self.log.getEffectiveLevel() <= logging.DEBUG and self.plot:
            for i,step in enumerate(steps):
                self.save(step,"%04d-Intermediate-%d: %s" % (lenslet,i,Stages[i]))
                self.frame().header=dict(Lenslet=lenslet,Spectrum=spectrum.label,Stage=Stages[i])
                self.frame().metadata=dict(Corner=corner)
                plt.imshow(step)
                plt.title("%s Image Generation Steps for Lenslet %4d" % (Stages[i],lenslet))
                plt.savefig("%s%04d-%s%s" % (self.config["System"]["Dirs"]["Partials"],lenslet,StagesF[i],self.config["Instrument"]["plot_format"]))
                plt.clf()
        
        # We only write the sub-image if the function is called to write sub images
        if write:
            self.write("%sSubimage-%4d%s" % (self.config["System"]["Dirs"]["Images"],lenslet,".fits"),clobber=True)
            self.keep(None)
        
        if do_return:
            return small, corner
        return
    
    
    # Placement Functions
    def place_cached_sed(self,lenslet,label,dlabel):
        """Places a cached SED Subimage"""
        slabel = "SUBIMG%d" % lenslet
        try:
            subframe = self.frame(slabel)
        except KeyError as e:
            raise SEDLimits(str(e))
        subimg = subframe()
        mlenslet = subframe.header['Lenslet']
        mcorner = subframe.metadata['Corner']
        if mlenslet != lenslet:
            raise ValueError("Lenslet Number Mismatch %d:%d for state %s in %s" % (lenslet,mlenslet,slabel,self))
        self.place(subimg,mcorner,label,dlabel)
    
    def place_sed(self,lenslet,spectrum,label,dlabel):
        """Place an SED based on a lenslet number and a spectrum"""
        small, corner = self.get_sub_image(lenslet,spectrum,fast=True)
        # Place uses the corner position to insert a sub-image into the master image
        self.place(small,corner,label,dlabel)
    
    def place(self,img,corner,label,dlabel):
        """Place the given AstroObject.AnalyticSpectra.AnalyticSpectrum onto the SEDMachine Image"""
        
        xstart = corner[0]
        xend = xstart + img.shape[0]
        ystart = corner[1]
        yend = ystart + img.shape[1]
        data = self.data(dlabel)
        
        if data.shape[0] < xend or data.shape[1] < yend:
            raise SEDLimits
        
        if xstart < 0 or ystart < 0:
            raise SEDLimits
        
        if xend < 0 or yend < 0:
            raise SEDLimits
        
        data[xstart:xend,ystart:yend] += img
        
        self.save(data,label)
        self.remove(dlabel)
        self.save(data,dlabel)
        self.log.debug("Placed a new image %s into %s at corner %s" % (label,dlabel,corner))
    
    
    # Basic Manipulation Functions
    def generate_blank(self):
        """Generates a blank SEDMachine Image"""
        self.save(np.zeros((self.config["Instrument"]["image_size"]["px"],self.config["Instrument"]["image_size"]["px"])).astype(np.int16),"Blank")
    
    def crop(self,x,y,xsize,ysize=None,label=None):
        """Crops the provided image to twice the specified size, centered around the x and y coordinates provided."""
        if not ysize:
            ysize = xsize
        cropped = self.states[self.statename].data[x-xsize:x+xsize,y-ysize:y+ysize]
        self.log.debug("Cropped and Saved Image")
        if label == None:
            label = "Cropped"
        self.remove(label)
        self.save(cropped,label)
    
    def generateGaussNoise(self,label=None,mean=10,std=2.0):
        """Generates a gaussian noise mask, saving to this object"""
        distribution = np.random.normal
        shape = (self.config["Instrument"]["ccd_size"]["px"],self.config["Instrument"]["ccd_size"]["px"])
        if label == None:
            label = "Gaussian Noise Mask (%2g,%2g)" % (mean,std)
        arguments = (mean,std,shape)
        noise = distribution(*arguments)
        
        self.save(noise,label)
    
    def generatePoissNoise(self,label=None,lam=2.0):
        """Generates a poisson noise mask, saving to this object"""
        distribution = np.random.poisson
        shape = (self.config["Instrument"]["ccd_size"]["px"],self.config["Instrument"]["ccd_size"]["px"])
        if label == None:
            label = "Poisson Noise Mask (%2g)" % (lam)
        arguments = (lam,shape)
        noise = distribution(*arguments)
        
        self.save(noise,label)
    
    def ccdCrop(self):
        """Crops the image to the appropriate ccd size"""
        x,y = self.center
        size = self.config["Instrument"]["ccd_size"]["px"] / 2.0
        self.crop(x,y,size,label=self.statename)
    
    def setupNoise(self):
        """Makes noise masks"""
        self.generatePoissNoise("Dark",self.config["Instrument"]["dark"]*self.config["Instrument"]["exposure"])
        self.generatePoissNoise("Bias",self.config["Instrument"]["bias"])
    
    def applyNoise(self,target):
        """Apply the noise masks to the target image label"""
        
        dark = self.data("Dark")
        bias = self.data("Bias")
        
        data = self.data(target)
        
        data += dark + bias
        
        self.remove(target)
        self.save(data,target)
    



