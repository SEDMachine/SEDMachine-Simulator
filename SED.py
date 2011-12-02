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

import AstroObject
from AstroObject.AstroSpectra import SpectraObject
from AstroObject.AstroImage import ImageObject,ImageFrame
from AstroObject.AnalyticSpectra import BlackBodySpectrum, AnalyticSpectrum, FlatSpectrum
from AstroObject.Utilities import *

from Utilities import *

__version__ = file("VERSION",'r').read()
__all__ = ["update","SEDLimits","Model"]

def update(d, u):
    """A deep update command for dictionaries.
    This is because the normal dictionary.update() command does not handle nested dictionaries."""
    for k, v in u.iteritems():
        if isinstance(v, collections.Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d


class SEDLimits(Exception):
    """A Basic Error-Differentiation Class.
    This error is used to express the fact that the SEDModel has encountered a spectrum which can't be placed as some part of it falls outside of the limits of the SED system."""
    pass


class Model(ImageObject):
    """This is a model container for the SEDMachine data simulator. This class is based on `AstroObject.ImageObject`, which it uses to provide some awareness of the way images are stored and retrieved. As such, this object has .save(), .select() and .show() etc. methods which can be used to examine the underlying data. It also means that this ImageObject subclass will contain the final, simulated image when the system is done.
    
    Note:: This object is not yet thread-safe, but coming versions of AstroObject will allow objects like this one to be locking and thread safe. Unfortunately, this limitation means that the simulator can generally not be run in multi-thread mode yet."""
    
    def __init__(self,configFile,scriptConfig):
        super(Model, self).__init__()
        self.configFile = configFile
        self.cache = False
        self.plot = True
        self.configs = {}
        self.config = {}
        if scriptConfig != None:
            self.scriptConfig = scriptConfig
            self.cache = scriptConfig["Cache"]
        self.regenerate = not self.cache
        self.WLS = {}
        self.initLog()
        self.configure()
    
    def initLog(self):
        """Setup Logging Functions for the SEDMachine Model.
        The logging system uses a Logs/ folder to store logs. If this folder is not present, it will not store logs.
        
        Once the logger is done initializing and configuring, it will print a starting message to the log file to differentiate between different simulation runs.
        
        Note:: There is a configuration file directive for the Log folder, but it is not respected at this point, as that would require pre-loading the configuration, without the ability to log that the configuration loading failed.
        """
        
        self.log = logging.getLogger(__name__)
        
        logfolder = "Logs/"
        filename = __name__+"-"+time.strftime("%Y-%m-%d")
        longFormat = "%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s"
        shortFormat = '%(levelname)-8s: %(message)s'
        dateFormat = "%Y-%m-%d-%H:%M:%S"
        
        self.log.setLevel(logging.DEBUG)
        logging.captureWarnings(True)
        self.console = logging.StreamHandler()
        self.console.setLevel(logging.INFO)
        self.consoleFormatter = logging.Formatter(shortFormat,datefmt=dateFormat)
        self.console.setFormatter(self.consoleFormatter)
        self.log.addHandler(self.console)
        
        if os.access(logfolder,os.F_OK):
            self.logfile = logging.FileHandler(filename=logfolder+filename+".log",mode="a")
            self.logfile.setLevel(logging.DEBUG)
            fileformatter = logging.Formatter(longFormat,datefmt=dateFormat)
            self.logfile.setFormatter(fileformatter)
            self.log.addHandler(self.logfile)
            self.log.removeHandler(self.console)
        
        self.log.info("--------------------------------")
        self.log.info("Welcome to the SED Machine model")
        self.log.info(" Version %s" % __version__ )
    
    
    # CONFIGURATION
    def configure(self):
        """This is a wrapper function for a variety of private configuration steps. The configuration generall happnes in stages:
        - Default Values are initalized
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
        self.config["convert"] = {}
        self.config["ccd_size"] = {}
        self.config["image_size"] = {}
        self.config["tel_radii"] = {}
        self.config["tel_obsc"] = {}
        self.config["psf_stdev"] = {}
        self.config["psf_size"] = {}
        # Unit Conversion Defaults
        self.config["convert"]["pxtomm"] = 0.0135
        # CCD / Image Plane Information
        self.config["ccd_size"]["px"] = 2048 #pixels
        self.config["image_size"]["mm"] = 40.0 #mm
        # Telescope Information
        self.config["tel_radii"]["px"] = 2.4 / 2.0
        self.config["tel_obsc"]["px"] = 0.4 / 2.0
        # PSF Information
        self.config["psf_size"]["px"] = 0
        # For a gaussian PSF
        self.config["psf_stdev"]["px"] = 1.0
        # Image Generation Density
        self.config["density"] = 5
        self.config["padding"] = 5
        # Default Gain Value
        self.config["gain"] = 1e-6
        # Noise Information
        self.config["dark"] = 20 # counts per pixel per second at some fixed degree c
        self.config["bias"] = 20 # counts per pixel at some fixed degree c
        self.config["exposure"] = 120 #Seconds
        # File Information for data and Caches
        self.config["files"] = {}
        self.config["files"]["lenslets"] = "Data/xy_17nov2011_v57.TXT"
        self.config["files"]["dispersion"] = "Data/dispersion_12-10-2011.txt"
        self.config["files"]["encircledenergy"] = "Data/encircled_energy_4nov11.TXT"
        # MPL Plotting Save Format
        self.config["plot_format"] = ".pdf"
        
        self.configs["defaults"] = copy.deepcopy(self.config)
    
    def _configureFile(self):
        """Attempt to load the default configuration from the working directory"""
        try:
            stream = file(self.configFile,'r')
        except IOError:
            self.log.warning("Configuration File Not Found: %s" % self.configFile)
        else:
            update(self.config,yaml.load(stream))
        finally:
            self.configs["file"] = copy.deepcopy(self.config)
    
    def _configureCaches(self):
        """If we are using the caching system, set up the configuration cache"""
        if not self.cache:
            self.log.debug("Skipping Caches")
            return
        
        self.log.debug("Enabling Caching")
        
        self.scriptConfig["Cache-Files"] = {}
        self.scriptConfig["Cache-Files"]["telescope"] = self.scriptConfig["Dirs"]["Caches"] + self.configFile.rstrip(".yaml") + ".tel"
        self.scriptConfig["Cache-Files"]["psf"] = self.scriptConfig["Dirs"]["Caches"] + self.configFile.rstrip(".yaml") + ".psf"
        self.scriptConfig["Cache-Files"]["conv"] = self.scriptConfig["Dirs"]["Caches"] + self.configFile.rstrip(".yaml") + ".conv"
        self.scriptConfig["Cache-Files"]["config"] = self.scriptConfig["Dirs"]["Caches"] + self.configFile
        self.scriptConfig["Cache-Files"]["wls"] = self.scriptConfig["Dirs"]["Caches"] + self.configFile.rstrip(".yaml") + ".wls"
        
        self.configs["precached"] = copy.deepcopy(self.config)
        
        try:
            stream = file(self.scriptConfig["Cache-Files"]["config"],'r')
            update(self.config,yaml.load(stream))
        except IOError:
            self.log.info("Cached Configuration File Not Found: %s" % self.scriptConfig["Cache-Files"]["config"])
        except TypeError:
            self.log.critical("Cached Configuration File has a Problem... skipping.")
        
        if self.configs["precached"] == self.config:
            self.regenearte = False
            self.log.debug("Configuration has not changed, will not regenerate kernels.")
        else:
            self.regenearte = True
            self.log.info("Configuration appears to have changed, will regenerate kernels.")
        
        
        return
    
    def _configureDynamic(self):
        """docstring for configureDynamicValues"""
        self.configs["NoDynamic"] = copy.deepcopy(self.config)
        
        if "calc" not in self.config["convert"]:
            if "mmtopx" not in self.config["convert"] and "pxtomm" in self.config["convert"]:
                self.config["convert"]["mmtopx"] = 1.0 / self.config["convert"]["pxtomm"]
                self.config["convert"]["calc"] = True
            else:
                self.config["convert"]["pxtomm"] = 1.0 / self.config["convert"]["mmtopx"]
                self.config["convert"]["calc"] = True
        
        self.config = self._setUnits(self.config,None)
        
        self.config["image_size"]["px"] = np.round( self.config["image_size"]["px"] , 0 )
    
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
                    r["px"] = v * self.config["convert"]["mmtopx"]
                    r["calc"] = True
                else:
                    self.log.warning("Value for %s set in both px and mm." % parent)
            elif k == "px":
                if ("calc" in r) and ("mm" in r):
                    pass
                elif ("calc" not in r):
                    r["mm"] = v * self.config["convert"]["pxtomm"]
                    r["calc"] = True
                else:
                    self.log.warning("Value for %s set in both px and mm." % parent)
        return r
    
    
    # Cacheing Functions
    def regenerateCache(self):
        """Cache calculated components of the system, including the telescope image and encircled energy image. Caches are stored to speed up system initalization. This function regenerates all cached files, including the configuration file. You can force the script to ignore cached files in the runner script using the option `--no-cache`. To regenrate the cache manually, simply delete the contents of the Caches directory"""
        stream = file(self.scriptConfig["Cache-Files"]["config"],'w')
        yaml.dump(self.configs["NoDynamic"],stream,default_flow_style=False)
        np.save(self.scriptConfig["Cache-Files"]["telescope"],self.TELIMG)
        np.save(self.scriptConfig["Cache-Files"]["psf"],self.PSFIMG)
        np.save(self.scriptConfig["Cache-Files"]["conv"],self.FINIMG)
    
    def cachedKernel(self):
        """Load cached kernels from the Caches directory. If any file is missing, it will attempt to trigger regeneration of the cache.
        You can force the script to ignore cached files in the runner script using the option `--no-cache`. To regenrate the cache manually, simply delete the contents of the Caches directory"""
        try:
            # Telescope Image Setup
            self.TELIMG = np.load(self.scriptConfig["Cache-Files"]["telescope"]+".npy")
            # PSF Setup
            self.PSFIMG = np.load(self.scriptConfig["Cache-Files"]["psf"]+".npy")
            # Preconvolved System
            self.FINIMG = np.load(self.scriptConfig["Cache-Files"]["conv"]+".npy")
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
            cachedWL = np.load(self.scriptConfig["Cache-Files"]["wls"]+".npy")
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
        stream = file(self.configFile,'w')
        yaml.dump(self.configs["NoDynamic"],stream,default_flow_style=False)
    
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
            self.log.debug("Generating Kernel Images")
            plt.clf()
            plt.imshow(self.TELIMG)
            plt.title("Telescope Image")
            plt.colorbar()
            plt.savefig("Partials/TelImage%s" % (self.config["plot_format"]))
            plt.clf()
            plt.imshow(self.PSFIMG)
            plt.title("PSF Image")
            plt.colorbar()
            plt.savefig("Partials/PSFImage%s" % (self.config["plot_format"]))
            plt.clf()
            plt.imshow(self.FINIMG)
            plt.title("Convolved Tel + PSF Image")
            plt.colorbar()
            plt.savefig("Partials/FINImage%s" % (self.config["plot_format"]))
            plt.clf()
        
        # Get layout data from files
        self.loadOpticsData(self.config["files"]["lenslets"],self.config["files"]["dispersion"])
        
        
        
        
        self.generate_blank()
        
        self.log.info("Done with SEDM setup")
    
    
    # Kernel Creation for Image Manipulation
    def regenerateKernel(self):
        """Regenerate a kernel from the configuration valuse"""
        self.log.info("Generating Kernels for Image System")
        # Telescope Image Setup
        self.TELIMG = self.get_tel_kern()
        # PSF Setup
        self.PSFIMG = self.get_psf_kern()
        # Preconvolved System
        self.FINIMG = sp.signal.convolve(self.PSFIMG,self.TELIMG,mode='same')
    
    def psf_kern(self,filename,size=0,truncate=False):
        """Generates a PSF Kernel from a file with mm-encircled energy conversions"""
        
        uM,FR = np.genfromtxt(filename,skip_header=18).T
        
        PX = uM * 1e-3 * self.config["convert"]["mmtopx"] * self.config["density"]
        
        if np.max(PX) <= size or truncate:
            size = np.int(size)
        else:
            size = np.int(np.max(PX))
        
        fit_vars = sp.interpolate.splrep(PX,FR)
        
        fit = lambda x : sp.interpolate.splev(x,fit_vars,der=1)
        
        vfit = np.vectorize(fit)
        
        x , y = np.mgrid[-size:size+1,-size:size+1]
        
        r = np.sqrt(x**2 + y**2)
        
        v = vfit(r)
        
        val = v
        
        self.log.debug("Generated a PSF Kernel for the encircled energy file %s with shape %s" % (filename,str(v.shape)))
        return val / np.sum(val)
    
    def circle_kern(self,radius,size=0,normalize=False):
        """Generate a Circle Kernel"""
        # This allows us to set the size of the kernel arbitrarily
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
        """ Returns a normalized 2D gauss kernel array for convolutions """
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
        """Returns the telescope kernel"""
        TELIMG = self.circle_kern( self.config["tel_radii"]["px"] * self.config["density"] )
        center = self.circle_kern( self.config["tel_obsc"]["px"] * self.config["density"] ,
            self.config["tel_radii"]["px"] * self.config["density"] , False )
        TELIMG -= center
        TELIMG = TELIMG / np.sum(TELIMG)
        self.log.debug("Generated a Telescpe Kernel with shape %s" % (str(TELIMG.shape)))
        return TELIMG
    
    def get_psf_kern(self):
        """Returns the PSF Kernel"""
        try:
            if self.config["psf_size"]["px"] != 0:
                size = self.config["psf_size"]["px"] * self.config["density"]
                truncate = True
            else:
                size = 0
                truncate = False
            
            PSFIMG = self.psf_kern(self.config["files"]["encircledenergy"],size,truncate)
        except IOError as e:
            self.log.warning("Could not access encircled energy file: %s" % e)
            PSFIMG = self.gauss_kern( (self.config["psf_stdev"]["px"] * self.config["density"]) )
        else:
            self.log.debug("Loaded Encircled Energy from %s" % self.config["files"]["encircledenergy"])
        return PSFIMG
    
    
    def get_blank_img(self):
        """Returns an image of the correct size to use for spectrum placement"""
        return np.zeros((self.config["image_size"]["px"],self.config["image_size"]["px"]))
    
    
    # Wavelength Functions
    def get_wavelengths(self,lenslet_num):
        """docstring for get_wavelengths"""
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
        """docstring for positionCaching"""
        if self.cache and self.regenerate:
            toCache = np.array(zip(self.WLS.keys(),self.WLS.values()))
            np.save(self.scriptConfig["Cache-Files"]["wls"],toCache)
    
    def _get_wavelengths(self,lenslet_num):
        """Returns an array of wavelengths for a given lenslet number"""
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
        npix = (distance * self.config["convert"]["mmtopx"]).astype(np.int) * self.config["density"]
        
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
        if self.config["density"] == 1:
            superDense_pts = superDense_pts.astype(np.int)
        else:
            superDense_pts = np.round(superDense_pts * self.config["density"]) / self.config["density"]
            superDense_int = (superDense_pts * self.config["density"]).astype(np.int)
        
        # We can identify unique points using the points when the integer position ratchets up or down.
        unique_x,unique_y = np.diff(superDense_int).astype(np.bool)
        
        # We want unique index to include points where either 'y' or 'x' ratchets up or down
        unique_idx = np.logical_or(unique_x,unique_y)
        
        # Remove any duplicate points. This does not do so in order, so we must
        # sort the array of unique points afterwards...
        unique_pts = superDense_pts[:,1:][:,unique_idx]
        
        # An array of distances to the origin of this spectrum, can be used to find wavelength
        # of light at each distance
        distance = superDense_distance[unique_idx] * self.config["convert"]["pxtomm"]
        
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
            plt.plot(points[:-1,1],np.diff(distance) * self.config["convert"]["mmtopx"],'g.')
            plt.savefig("Partials/%04d-Distances_Diff%s" % (lenslet_num,self.config["plot_format"]))
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
        """docstring for Load Lenselt Data"""
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
        xs += (self.config["image_size"]["mm"]/2)
        ys += (self.config["image_size"]["mm"]/2)
        
        # Find the xs and ys that are not within 0.1 mm of the edge of the detector...
        ok = (xs > 0.1) & (xs < self.config["image_size"]["mm"]-0.1) & (ys > 0.1) & (ys < self.config["image_size"]["mm"]-0.1)
        ix, p1, p2, lams, xs, ys = ix[ok], p1[ok], p2[ok], lams[ok], xs[ok], ys[ok]
        # We remove these positions because we don't want to generate spectra for them.
        
        # This simply generates a list of all of the lenslets
        self.lenslets = np.unique(ix)
        
        # Convert the xs and ys to pixel positions
        xpix = np.round(xs * self.config["convert"]["mmtopx"],0).astype(np.int)
        ypix = np.round(ys * self.config["convert"]["mmtopx"],0).astype(np.int)
        
        # Determine the center of the whole system by finding the x position that is closest to 0,0 in pupil position
        cntix = np.argmin(p1**2 + p2**2)
        self.center = (xs[cntix] * self.config["convert"]["mmtopx"], ys[cntix] * self.config["convert"]["mmtopx"])
        
        self.ix, self.p1, self.p2, self.lams, self.xs, self.ys = ix, p1, p2, lams, xs, ys
        self.xpix, self.ypix = xpix, ypix
    
    def loadOpticsData(self,laspec,dispspec):
        """Loads an optical conversion based on the lenslet array spec and dispersion spec files provided"""
        
        self.loadDispersionData(dispspec)
        self.loadLensletData(laspec)
    
    
    # Image Tracing Function
    def get_dense_image(self,lenslet,spectrum):
        """This function returns a dense image array, with flux placed into single pixels."""
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
        radiance = spectrum(wl[:-1]*1e-6) * self.config["gain"]
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
            plt.savefig("Partials/%04d-SpecFlux%s" % (lenslet,self.config["plot_format"]))
            plt.clf()
            plt.plot(wl[:-1],deltawl,"g.")
            plt.title("$\Delta\lambda$ for each pixel")
            plt.xlabel("Wavelength ($\mu m$)")
            plt.ylabel("$\Delta\lambda$ per pixel")
            plt.savefig("Partials/%04d-SpecDeltaWL%s" % (lenslet,self.config["plot_format"]))
            plt.clf()
        
        # Take our points out. Note from the above that we multiply by the density in order to do this
        xorig,yorig = (points * self.config["density"])[:-1].T.astype(np.int)
        x,y = (points * self.config["density"])[:-1].T.astype(np.int)
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
        xdist += (self.config["density"] - xdist % self.config["density"])
        ydist += (self.config["density"] - ydist % self.config["density"])
        
        # Move our x and y coordinates to the middle of our sub image by applying padding below each one.
        x += self.config["padding"] * self.config["density"]
        y += self.config["padding"] * self.config["density"]
        
        # Find the first (by the flatten method) corner of the subimage,
        # useful for placing the sub-image into the full image.
        corner = np.array([ xint[np.argmax(x)], yint[np.argmin(y)]])
        self.log.debug("Corner Position in Integer Space: %s" % corner)
        corner *= self.config["density"]
        realcorner = np.array([ xorig[np.argmax(x)], yorig[np.argmin(y)]])
        offset = corner - realcorner
        corner /= self.config["density"]
        self.log.debug("Corner Position Offset in Dense Space: %s" % (offset))
        if self.log.getEffectiveLevel() <= logging.DEBUG:
            handle = file("Partials/Offsets.dat",'a')
            np.savetxt(handle,offset.T)
        corner -= np.array([-self.config["padding"],self.config["padding"]])
        
        x += offset[0]
        y += offset[1]
        
        # Create our sub-image, using the x and y width of the spectrum, plus 2 padding widths.
        # Padding is specified in full-size pixels to ensure that the final image is an integer
        # number of full-size pixels across.
        img = np.zeros((xdist+2*self.config["padding"]*self.config["density"],ydist+2*self.config["padding"]*self.config["density"]))
        
        # Place the spectrum into the sub-image
        img[x,y] = flux
        self.log.debug("Placing flux (shape: %s ) for spectrum %4d into a sub-image (shape: %s)."
            % (flux.shape,lenslet,img.shape))
        
        if self.log.getEffectiveLevel() <= logging.DEBUG:
            np.savetxt("Partials/SubimageValues.dat",np.array([x,y,wl[:-1],deltawl,flux]).T)
        
        return img, corner
    
    
    # Image Convolution
    def get_sub_image(self,lenslet,spectrum,fast=False):
        """Returns a sub-image for a given lenslet"""
        # This function gets the dense sub-image with the spectrum placed in single dense pixels
        img,corner = self.get_dense_image(lenslet,spectrum)
        self.log.debug("Retrieved Dense Image for %4d" % lenslet)
        
        if fast:
            # Convolve with the PSF and Telescope Image simultaneously
            img2 = sp.signal.convolve(img,self.FINIMG,mode='same')
            self.log.debug("Convolved Dense Image with PSF and Telescope for %4d" % lenslet)
            # Bin the image back down to the final pixel size
            small = bin(img2,self.config["density"]).astype(np.int16)
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
            small = bin(img_tel_psf,self.config["density"]).astype(np.int16)
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
                plt.savefig("Partials/%04d-%s%s" % (lenslet,StagesF[i],self.config["plot_format"]))
                plt.clf()
        
        # We only write the sub-image if the function is called to write sub images
        if write:
            self.write("Images/Subimage-%4d%s" % (lenslet,".fits"),clobber=True)
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
    
    
    # Basic Manipulation Functions
    def generate_blank(self):
        """Generates a blank SEDMachine Image"""
        self.save(np.zeros((self.config["image_size"]["px"],self.config["image_size"]["px"])).astype(np.int16),"Blank")
    
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
        shape = (self.config["ccd_size"]["px"],self.config["ccd_size"]["px"])
        if label == None:
            label = "Gaussian Noise Mask (%2g,%2g)" % (mean,std)
        arguments = (mean,std,shape)
        noise = distribution(*arguments)
        
        self.save(noise,label)
    
    def generatePoissNoise(self,label=None,lam=2.0):
        """Generates a poisson noise mask, saving to this object"""
        distribution = np.random.poisson
        shape = (self.config["ccd_size"]["px"],self.config["ccd_size"]["px"])
        if label == None:
            label = "Poisson Noise Mask (%2g)" % (lam)
        arguments = (lam,shape)
        noise = distribution(*arguments)
        
        self.save(noise,label)
    
    def ccdCrop(self):
        """Crops the image to the appropriate ccd size"""
        x,y = self.center
        size = self.config["ccd_size"]["px"] / 2.0
        self.crop(x,y,size,label=self.statename)
    
    def setupNoise(self):
        """Makes noise masks"""
        self.generatePoissNoise("Dark",self.config["dark"]*self.config["exposure"])
        self.generatePoissNoise("Bias",self.config["bias"])
    
    def applyNoise(self,target):
        """Apply the noise masks to the target image label"""
        
        dark = self.data("Dark")
        bias = self.data("Bias")
        
        data = self.data(target)
        
        data += dark + bias
        
        self.remove(target)
        self.save(data,target)
    



