#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SED.py
#  Simulation Software
#
#  Created by Alexander Rudy on 2011-11-04.
#  Copyright 2011 Alexander Rudy. All rights reserved.
#  Version 0.2.0a1
#


import numpy as np
import pyfits as pf
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.signal
import scipy.interpolate
import yaml

import arpytools.progressbar

import os,logging,time,copy,collections

import logging.handlers

import AstroObject
import AstroObject.AstroSimulator
from AstroObject.AstroSpectra import SpectraObject
from AstroObject.AstroImage import ImageObject,ImageFrame
from AstroObject.AnalyticSpectra import BlackBodySpectrum, AnalyticSpectrum, FlatSpectrum
from AstroObject.Utilities import *

from Lenslet import *
from Utilities import *

__version__ = getVersion(__name__)
__all__ = ["SEDLimits","Instrument"]
LOG = logging.getLogger(__name__)


class SubImage(ImageFrame):
    """A custom frame type for SEDMachine Sub-Images"""
    def __init__(self, array, label, header=None, metadata=None):
        super(SubImage, self).__init__(array,label,header,metadata)
        self.lensletNumber = 0
        self.corner = [0,0]
        self.configHash = hash(0)
        self.spectrum = "NO SPEC"
        
    def sync_header(self):
        """Synchronizes the header dictionary with the HDU header"""
        # assert self.label == self.header['SEDlabel'], "Label does not match value specified in header: %s and %s" % (self.label,self.header['SEDlabel'])
        self.configHash = self.header['SEDconf']
        self.corner = [self.header['SEDcrx'],self.header['SEDcry']]
        self.spectrum = self.header['SEDspec']
        self.lensletNumber = self.header['SEDlens']
    
    def __hdu__(self,primary=False):
        """Retruns an HDU which represents this frame. HDUs are either ``pyfits.PrimaryHDU`` or ``pyfits.ImageHDU`` depending on the *primary* keyword."""
        if primary:
            LOG.log(5,"Generating a primary HDU for %s" % self)
            HDU = pf.PrimaryHDU(self())
        else:
            LOG.log(5,"Generating an image HDU for %s" % self)
            HDU = pf.ImageHDU(self())
        HDU.header.update('object',self.label)
        HDU.header.update('SEDlabel',self.label)
        HDU.header.update('SEDconf',self.configHash)
        HDU.header.update('SEDcrx',self.corner[0])
        HDU.header.update('SEDcry',self.corner[1])
        HDU.header.update('SEDspec',self.spectrum)
        HDU.header.update('SEDlens',self.lensletNumber)
        for key,value in self.header.iteritems():
            HDU.header.update(key,value)
        return HDU
    
    def __show__(self):
        """Plots the image in this frame using matplotlib's ``imshow`` function. The color map is set to an inverted binary, as is often useful when looking at astronomical images. The figure object is returned, and can be manipulated further.
        
        .. Note::
            This function serves as a quick view of the current state of the frame. It is not intended for robust plotting support, as that can be easily accomplished using ``matplotlib``. Rather, it attempts to do the minimum possible to create an acceptable image for immediate inspection.
        """
        LOG.debug("Plotting %s using matplotlib.pyplot.imshow" % self)
        figure = plt.imshow(self())
        figure.set_cmap('binary_r')
        return figure
    
    @classmethod
    def __read__(cls,HDU,label):
        """Attempts to convert a given HDU into an object of type :class:`ImageFrame`. This method is similar to the :meth:`__save__` method, but instead of taking data as input, it takes a full HDU. The use of a full HDU allows this method to check for the correct type of HDU, and to gather header information from the HDU. When reading data from a FITS file, this is the prefered method to initialize a new frame.
        """
        LOG.debug("Attempting to read as %s" % cls)
        if not isinstance(HDU,(pf.ImageHDU,pf.PrimaryHDU)):
            msg = "Must save a PrimaryHDU or ImageHDU to a %s, found %s" % (cls.__name__,type(HDU))
            raise AbstractError(msg)
        if not isinstance(HDU.data,np.ndarray):
            msg = "HDU Data must be %s for %s, found data of %s" % (np.ndarray,cls.__name__,type(HDU.data))
            raise AbstractError(msg)
        try:
            Object = cls(HDU.data,label,header=HDU.header)
        except AssertionError as AE:
            msg = "%s data did not validate: %s" % (cls.__name__,AE)
            raise AbstractError(msg)
        LOG.debug("Created %s" % Object)
        Object.sync_header()
        return Object
    




class Instrument(ImageObject,AstroObject.AstroSimulator.Simulator):
    """This is a model container for the SEDMachine data simulator. This class is based on `AstroObject.ImageObject`, which it uses to provide some awareness of the way images are stored and retrieved. As such, this object has .save(), .select() and .show() etc. methods which can be used to examine the underlying data. It also means that this ImageObject subclass will contain the final, simulated image when the system is done.
    
    .. Note:: This object is not yet thread-safe, but coming versions of AstroObject will allow objects like this one to be locking and thread safe. Unfortunately, this limitation means that the simulator can generally not be run in multi-thread mode yet.
    
    """
    
    def __init__(self,config):
        super(Instrument, self).__init__(name="Instrument")
        self.dataClasses = [SubImage]
        self.config = config
        self.debug = self.config["Debug"]
        self.cache = self.config["Cache"]
        self.plot = self.config["Plot"]
        self.defaults = []
        self.regenerate = not self.cache
        self.WLS = {}
        self._configure()
    
    def update(self, d, u):
        """A deep update command for dictionaries.
        This is because the normal dictionary.update() command does not handle nested dictionaries.
        """
        for k, v in u.iteritems():
            if isinstance(v, collections.Mapping):
                r = self.update(d.get(k, {}), v)
                d[k] = r
            else:
                d[k] = u[k]
        return d
    
    
    
    # CONFIGURATION
    def _configure(self):
        """This is a wrapper function for a variety of private configuration steps. The configuration generall happnes in stages:
        
        - Default Values are initalized for the instrument configuration
        - System loads the configuration file, updating the default values with the user-selected values. If no configuration file exists, the system will move on.
        - System loads the configuration from the Caches directory if caching is enabled.
        - System adds dynamic variables to the configuration. This is the feature which handles variable units etc.
        
        To see what the default configuration file looks like, run the runner script with the --dump argument. The runner script will requrie a subcommand, but dump will cause the program to exit before recieving any subcommands for action.
        
        You can force the script to ignore cached files in the runner script using the option `--no-cache`. To regenrate the cache manually, simply delete the contents of the Caches directory"""
        self._configureDefaults()
        self.config["Configs"]["This"] = self.config["Configs"]["Instrument"]
        self.configure(configuration=self.config)
        self._configureCaches()
        self._configureDynamic()
        fileName = "%(dir)sInstrument-Config.dat" % {'dir':self.config["Dirs"]["Partials"]}
        with open(fileName,'w') as stream:
            yaml.dump(self.config,stream,default_flow_style=False)
    
    def _configureDefaults(self):
        """Set up the default configure variable. If you change the default configuration variables in this function (instead of using a configuration file), the script will generally not detect the change, and so will not regenerate Cached files. You can force the script to ignore cached files in the runner script using the option `--no-cache`. To regenrate the cache manually, simply delete the contents of the Caches directory
        """
        self.defaultConfig = {}
        # Configuration Variables for The System
        # Unit Conversion Defaults
        self.defaultConfig["convert"] = {}
        self.defaultConfig["convert"]["pxtomm"] = 0.0135
        # CCD / Image Plane Information
        self.defaultConfig["ccd_size"] = {}
        self.defaultConfig["ccd_size"]["px"] = 2048 #pixels
        self.defaultConfig["image_size"] = {}
        self.defaultConfig["image_size"]["mm"] = 40.0 #mm
        # Telescope Information
        self.defaultConfig["tel_radii"] = {}
        self.defaultConfig["tel_radii"]["px"] = 2.4 / 2.0
        self.defaultConfig["tel_obsc"] = {}
        self.defaultConfig["tel_obsc"]["px"] = 0.4 / 2.0
        # PSF Information
        self.defaultConfig["psf_size"] = {}
        self.defaultConfig["psf_size"]["px"] = 0
        # For a gaussian PSF
        self.defaultConfig["psf_stdev"] = {}
        self.defaultConfig["psf_stdev"]["px"] = 1.0
        # Image Generation Density
        self.defaultConfig["density"] = 5
        self.defaultConfig["padding"] = 5
        # Default Gain Value
        self.defaultConfig["gain"] = 1e10
        # Noise Information
        self.defaultConfig["dark"] = 20 # counts per pixel per second at some fixed degree c
        self.defaultConfig["bias"] = 20 # counts per pixel at some fixed degree c
        self.defaultConfig["exposure"] = 120 #Seconds
        # File Information for data and Caches
        self.defaultConfig["files"] = {}
        self.defaultConfig["files"]["lenslets"] = "Data/xy_17nov2011_v57.TXT"
        self.defaultConfig["files"]["dispersion"] = "Data/dispersion_12-10-2011.txt"
        self.defaultConfig["files"]["encircledenergy"] = "Data/encircled_energy_4nov11.TXT"
        # MPL Plotting Save Format
        self.defaultConfig["plot_format"] = ".pdf"
        self.defaultConfig["plot"] = False
        
        # Logging Configuration
        self.defaultConfig["logging"] = {}
        self.defaultConfig["logging"]["console"] = {}
        self.defaultConfig["logging"]["console"]["enable"] = True
        self.defaultConfig["logging"]["console"]["format"] = "......%(message)s"
        self.defaultConfig["logging"]["console"]["level"] = False
        self.defaultConfig["logging"]["file"] = {}
        self.defaultConfig["logging"]["file"]["enable"] = True
        self.defaultConfig["logging"]["file"]["filename"] = "SEDInstrument"
        self.defaultConfig["logging"]["file"]["format"] = "%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s"
        
        self.defaults += [copy.deepcopy(self.defaultConfig)]
        
        self.update(self.config,self.defaultConfig)
        
        self.log.debug("Set Instrument Configuration Defaults")
    
    def _configureCaches(self):
        """If we are using the caching system, set up the configuration cache.
        
        First, we set default Caching variables. Next, load these same variables (and overwrite) from the system configuration. Next, we load the cached configuration. If the cached configuration updates the generated configuration in any way, """
        if not self.cache:
            self.log.debug("Skipping Cache Configuration")
            return
        
        # Default Cache Variables
        CacheFiles = {}
        CacheFileBase = "SED.Instrument"
        CacheFiles["telescope"] = CacheFileBase + ".tel"    + ".npy"
        CacheFiles["psf"]       = CacheFileBase + ".psf"    + ".npy"
        CacheFiles["conv"]      = CacheFileBase + ".conv"   + ".npy"
        CacheFiles["config"]    = CacheFileBase + ".config" + ".yaml"
        CacheFiles["wls"]       = CacheFileBase + ".wls"    + ".npz"
        
        CacheFiles = self.update(CacheFiles,self.config["CacheFiles"])
        self.config["CacheFiles"] = CacheFiles
        
        self.log.debug("Configured Cache Variables")
        self.defaults += [copy.deepcopy(self.config)]
        self.Caches.directory = self.config["Dirs"]["Caches"]
        self.Caches.registerNPY("TEL",self.get_tel_kern,filename=self.config["CacheFiles"]["telescope"])
        self.Caches.registerNPY("PSF",self.get_psf_kern,filename=self.config["CacheFiles"]["psf"])
        self.Caches.registerNPY("CONV",lambda : sp.signal.convolve(self.Caches.get("PSF"),self.Caches.get("TEL"),mode='same'),filename=self.config["CacheFiles"]["conv"])
        self.Caches.registerCustom("CONFIG",kind=AstroObject.AstroSimulator.YAMLCache,generate=lambda : self.config,filename=self.config["CacheFiles"]["config"])
        
        
        if self.Caches.check(master="CONFIG"):
            self.log.info("Caches appear out of date, regenerating")
        self.Caches.get("CONFIG")
        
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
    
    def cachedWL(self):
        """Load cached wavelengths from the Caches directory. If any file is missing, it will attempt to trigger regeneration of the cache.
        You can force the script to ignore cached files in the runner script using the option `--no-cache`. To regenrate the cache manually, simply delete the contents of the Caches directory
        """
        try:
            # Cached Wavelengths
            self.WLS = {}
            ZIPFile = np.load(self.config["CacheFiles"]["wls"])
            for label in ZIPFile:
                self.WLS[label] = ZIPFile[label]
        except IOError as e:
            self.log.warning("Cached files not found, regenerating. Error: %s" % e )
            self.regenerate = True
        else:
            self.log.debug("Loaded Wavelengths from Numpy Files")
            with open("%s%s%s" % (self.config["Dirs"]["Partials"],"Instrument-CachedWls-Keys-Loaded",".dat"),'w') as stream:
                stream.write("\n".join(self.WLS.keys()))
    
    def resetWLCache(self):
        """Resets the Wavelength Cache, which will force it to regenerate during the simulation
        """
        self.log.debug("Forcing Wavelenghts to Regenerate")
        self.WLS = {}
    
    def dumpConfig(self):
        """Dumps a valid configuration file on top of any old configuration files. This is useful for examining the default configuration fo the system, and providing modifications through this interface."""
        with open(self.config["Configs"]["Instrument"].rstrip(".yaml")+".dump.yaml",'w') as stream:
            yaml.dump(self.defaults[1],stream,default_flow_style=False)
    
    def setup(self):
        """After the object has been initialized, it must be setup. Setup relies on an established configuration to determine what fixed parts of the system should be generated. Actions taken in the setup phase are:
        
        - Load data from the Cache (if enabled)
        - Regenerate Kernel Functions (if required)
        - Regenerate Cache Files (if required, and if enabled)
        - Load optical system data (not cached right now)
        - Generate a blank image (no need to cache... np.zeros is fast!)
        
        """
        
        if self.log.getEffectiveLevel() <= logging.DEBUG and self.config["plot"]:
            self.plotKernalPartials()
        
        # Get layout data from files
        self.loadOpticsData(self.config["files"]["lenslets"],self.config["files"]["dispersion"])
        self.setupLenslets()
        # TODO: We could cache this...
        
        self.generate_blank()
        
        self.log.info("Done with SEDM setup")
    
    
    def plotKernalPartials(self):
        """Plots the kernel data partials"""
        self.log.debug("Generating Kernel Plots and Images")
        plt.clf()
        plt.imshow(self.Caches.get("TEL"))
        plt.title("Telescope Image")
        plt.colorbar()
        plt.savefig("%sInstrument-TEL-Kernel%s" % (self.config["Dirs"]["Partials"],self.config["plot_format"]))
        plt.clf()
        plt.imshow(self.Caches.get("PSF"))
        plt.title("PSF Image")
        plt.colorbar()
        plt.savefig("%sInstrument-PSF-Kernel%s" % (self.config["Dirs"]["Partials"],self.config["plot_format"]))
        plt.clf()
        plt.imshow(self.Caches.get("CONV"))
        plt.title("Convolved Tel + PSF Image")
        plt.colorbar()
        plt.savefig("%sInstrument-FIN-Kernel%s" % (self.config["Dirs"]["Partials"],self.config["plot_format"]))
        plt.clf()
    
    def psf_kern(self,filename,size=0,truncate=False,header_lines=18):
        """Generates a PSF Kernel from a file with micron-encircled energy conversions. The file should have two columns, first, microns from the center of the PSF, and second, the fraction of encircled energy at that distance from the PSF.
        
        The calculation is then carried out using a spline fit to this data. From the spline fit, the function returns the first derivative of the encircled energy at each point. This in effect is the amount of energy at each point. These values are then normalized, to create a PSF mask for the instrument.
        
        The `size` parameter specifies the size of the kernel to use. If the size is greater than the encircled energy data, then a larger figure will be returned. If `size` is smaller than the encircled energy data, it will return an image the size of the encircled energy data, unless the `truncate` parameter is set to `true`.
        
        The `header_lines` parameter defaults to 18, which works with Zemax encircled energy output."""
        
        uM,FR = np.genfromtxt(filename,skip_header=header_lines).T
        # Convert microns to milimeters, then pixels, then dense pixels
        PX = uM * 1e-3 * self.config["convert"]["mmtopx"] * self.config["density"]
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
        TELIMG = self.circle_kern( self.config["tel_radii"]["px"] * self.config["density"] )
        center = self.circle_kern( self.config["tel_obsc"]["px"] * self.config["density"] ,
            self.config["tel_radii"]["px"] * self.config["density"] , False )
        TELIMG -= center
        TELIMG = TELIMG / np.sum(TELIMG)
        self.log.debug("Generated a Telescpe Kernel with shape %s" % (str(TELIMG.shape)))
        return TELIMG
    
    def get_psf_kern(self):
        """Returns the PSF Kernel. The function first tries to read the encircled energy file. In this case, if a `psf_size` is set in the instrument configuration, this value will be used to truncate the size of the encircled energy function. If the encircled energy function cannot be loaded, the system will fall back on to a gaussian psf as configured by the instrument."""
        if self.config["psf_size"]["px"] != 0:
            size = self.config["psf_size"]["px"] * self.config["density"]
            truncate = True
        else:
            size = 0
            truncate = False
        try:
            PSFIMG = self.psf_kern( self.config["files"]["encircledenergy"],size,truncate)
        except IOError as e:
            self.log.warning("Could not access encircled energy file: %s" % e)
            PSFIMG = self.gauss_kern( (self.config["psf_stdev"]["px"] * self.config["density"]) )
        else:
            self.log.debug("Loaded Encircled Energy from %s" % self.config["files"]["encircledenergy"])
        return PSFIMG
    
    
    def get_psf(self,wl):
        """Return a PSF for a given wavelength"""
        if wl > 5 * 1e-7:
            return self.Caches.get("CONV")
        else:
            return self.gauss_kern( (self.config["psf_stdev"]["px"] * self.config["density"] / 2) )
        
    def get_blank_img(self):
        """Returns an image of the correct size to use for spectrum placement. Set by the `image_size` value."""
        return np.zeros((self.config["image_size"]["px"],self.config["image_size"]["px"]))
    
    
    # Wavelength Functions
    def lenslet_dispersion(self):
        """Calculate the wavelengths for each lenslet"""
        
        PBar = arpytools.progressbar.ProgressBar(color="green")
        finished = 0.0
        total = len(self.lenslets)
        
        self.log.info("Calculating Wavelength Dispersion for %d lenslets" % total)
        
        self.log.useConsole(False)
        PBar.render(0,"L:%4d %4d/%-4d" % (0,finished,total))
        
        for index in self.lenslets:
            self.lensletObjects[index].find_dispersion()
            finished += 1.0
            progress = int((finished/float(total)) * 100)
            PBar.render(progress,"L:%4d %4d/%-4d" % (index,finished,total))
        
        self.log.useConsole(True)
        
        
    def lenslet_trace(self,spectrum):
        """Get trace for all lenslets"""
        PBar = arpytools.progressbar.ProgressBar(color="blue")
        finished = 0.0
        total = len(self.lenslets)
        
        self.log.info("Calculating Wavelength Trace for %d lenslets" % total)
        
        self.log.useConsole(False)
        PBar.render(0,"L:%4d %4d/%-4d" % (0,finished,total))
        
        for index in self.lenslets:
            self.lensletObjects[index].get_trace(spectrum(index))
            finished += 1.0
            progress = int((finished/float(total)) * 100)
            PBar.render(progress,"L:%4d %4d/%-4d" % (index,finished,total))
        
        self.log.useConsole(True)
        
          
            
        

    def lenslet_image(self):
        """Creates an image for each lenslet in desnse space"""
        PBar = arpytools.progressbar.ProgressBar(color="blue")
        finished = 0.0
        total = len(self.lenslets)
        
        self.log.info("Placing Wavelength Dependent PSFs for %d lenslets" % total)
        
        self.log.useConsole(False)
        PBar.render(0,"L:%4d %4d/%4d" % (0,finished,total))
        
        for index in self.lenslets:
            self.lensletObjects[index].place_trace(self.get_psf)
            finished += 1.0
            progress = int((finished/float(total)) * 100)
            PBar.render(progress,"L:%4d %4d/%-4d" % (index,finished,total))
        
        self.log.useConsole(True)
      
    def gen_hexagons(self):
        """Plot the hexagons"""
        PBar = arpytools.progressbar.ProgressBar(color="green")
        finished = 0.0
        total = len(self.lenslets)
        
        self.log.info("Generating Hexagons for %d lenslets" % total)
        
        self.log.useConsole(False)
        PBar.render(0,"L:%4d %4d/%4d" % (0,finished,total))
        
        for index in self.lenslets:
            self.lensletObjects[index].make_hexagon()
            finished += 1.0
            progress = int((finished/float(total)) * 100)
            PBar.render(progress,"L:%4d %4d/%-4d" % (index,finished,total))
        self.log.useConsole(True)
    
       
    def show_hexagons(self):
        """Plot the hexagons"""
        PBar = arpytools.progressbar.ProgressBar(color="cyan")
        finished = 0.0
        total = len(self.lenslets)
        
        self.log.info("Plotting Hexagons for %d lenslets" % total)
        
        
        self.log.useConsole(False)
        PBar.render(0,"L:%4d %4d/%4d" % (0,finished,total))
        plt.figure()
        plt.title("Hexagons")
        
        for index in self.lenslets:
            x, y = self.lensletObjects[index].shape.exterior.xy
            plt.fill(x, y, color='#cccccc', aa=True) 
            plt.plot(x, y, color='#666666', aa=True, lw=0.25)
            finished += 1.0
            progress = int((finished/float(total)) * 100)
            PBar.render(progress,"L:%4d %4d/%-4d" % (index,finished,total))
        plt.savefig("%sLenslets-Hexagon%s" % (self.config["Dirs"]["Partials"],self.config["plot_format"]))
        self.log.useConsole(True)
        
    
    def get_wavelengths(self,lenslet_num):
        """Returns a tuple of ([x,y],wavelength,delta wavelength). Thus, this function calculates the x,y pixel positions of the spectra, the wavelength at each pixel, and the delta wavelength for each pixel.
        
        .. Note:: This function can use cached calculations, which do not seem to necessarily speed up calculation.
        
        .. Note:: The logic of this collection is in `_get_wavelengths()`, which is heavily documented in source code.
        """
        
        
        if str(lenslet_num) in self.WLS:
            self.log.debug("Using Cached WLs for %d" % lenslet_num)
            xs,ys,wls = self.WLS[str(lenslet_num)]
            results = np.array([xs,ys]).T, wls, np.diff(wls)
        else:
            self.log.debug("Regenerating WLs for %d" % lenslet_num)
            results = self._get_wavelengths(lenslet_num)
            points,wls,wldiffs = results
            xs,ys = points.T
            self.WLS[str(lenslet_num)] = np.array([xs,ys,wls])
            if self.cache:
                for item,label in zip(results,["Points,WLs,DeltaWLs"]):
                    self.log.debug(npArrayInfo(item,label))
                self.regenerate |= True
        if len(results) != 3:
            self.log.critical("We have a problem with retrieved WLs")
            self.log.info(npArrayInfo(results,"WL Retrieved Results"))
            raise Exception
        
        filename = "%(dir)s%(file)s%(fmt)s" % { 'dir' : self.config["Dirs"]["Partials"],
            'file' : 'Instrument-Cached-WL-Points', 'fmt' : '.dat'}
        with open(filename,'w') as stream:
            stream.write(npArrayInfo(results[0],"Points") + "\n")
            for row in results[0]:
                for col in row:
                    stream.write(str(col)+" -> "+str(type(col))+" ")
                stream.write("| "+str(type(row))+"\n")
        self.log.debug(npArrayInfo(results[1],"CachedWLs"))
        return results
    
    def positionCaching(self):
        """Save the wavelength and position cache. This must be done after all wavelengths have been collected."""
        if self.cache and self.regenerate:
            np.savez(self.config["CacheFiles"]["wls"],**self.WLS)
        with open("%s%s%s" % (self.config["Dirs"]["Partials"],"Instrument-CachedWls-Keys-Saved",".dat"),'w') as stream:
            stream.write("\n".join(self.WLS.keys()))
        
    
    def _get_wavelengths(self,lenslet_num):
        """Returns an array of wavelengths for a given lenslet number. See `get_wavelenghts()`.
        
        .. Note:: This function is heavily documented in the source code.
        """
        # This method will be the critical slow point for the whole system...
        # It should be re-written to do the conversions in a better way, but
        # I'm not sure how to do that right now.
        
        # First, we take only data points which apply to this lenslet
        lenslet = self.lensletObjects[lenslet_num]
        # use = lenslet_num == self.ix
        
        # Interpolation to convert from wavelength to pixels.
        #   The accuracy of this interpolation is not important.
        #   Rather, it is used to find the pixels where the light will fall
        #   and is fed an array that is very dense, used on this dense interpolation
        #   and then binned back onto pixels. Thus it will be used to get a list
        #   of all illuminated pixels.
        fx = np.poly1d(np.polyfit(lenslet.ls, lenslet.xpixs, 2))
        fy = np.poly1d(np.polyfit(lenslet.ls, lenslet.ypixs, 2))
        
        # Find the starting and ending position of the spectra
        startix = np.argmin(lenslet.ls)
        endix = np.argmax(lenslet.ls)
        start = np.array([lenslet.xs[startix],lenslet.ys[startix]])
        end = np.array([lenslet.xs[endix],lenslet.ys[endix]])
        
        # Get the total length of the spectra
        distance = np.sqrt(np.sum(end-start)**2)
        
        if distance == 0:
            raise SEDLimits
        
        # Find the length in units of (int) pixels
        npix = (distance * self.config["convert"]["mmtopx"]).astype(np.int) * self.config["density"]
        
        # Create a data array one hundred times as dense as the number of pixels
        #   This is the super dense array which will use the above interpolation
        superDense_lam = np.linspace(np.min(lenslet.ls),np.max(lenslet.ls),npix*100)
        
        # Interpolate along our really dense set of wavelengths to find all possible
        # illuminated pixel positions in this spectrum
        superDense_pts = np.array([fx(superDense_lam),fy(superDense_lam)])
        
        # Measure the distance along our really dense set of points
        superDense_interval = np.sqrt(np.sum(np.power(np.diff(superDense_pts,axis=1),2),axis=0))
        superDense_distance = np.cumsum(superDense_interval)
        
        # Adjust the density of our points. This rounds all values to only full pixel values.
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
        self.log.debug(npArrayInfo(points,"Points"))
        
        
        # Pull out the original wavelengths
        wl_orig = superDense_lam[unique_idx][sorted_idx]
        wl = wl_orig
        self.log.debug(npArrayInfo(wl,"Wavelengths"))
        # We are getting some odd behavior, where the dispersion function seems to not cover the whole
        # arc length and instead covers only part of it. This causes much of our arc to leave the desired
        # and available wavelength boundaries. As such, I'm disabling the more accurate dispersion mode.
        
        # Convert to wavelength space along the dispersion spline.
        # wl = self.spline(distance)
        
        # This is a debugging area which will make a lot of graphs and annoying messages.
        if self.log.getEffectiveLevel() <= logging.DEBUG and self.config["plot"]:
            self.log.debug("Plotting Wavelength Information")
            # This graph shows the change in distance along arc per pixel.
            # The graph should produce all data points close to each other, except a variety of much lower
            # data points which are caused by the arc crossing between pixels.
            plt.clf()
            plt.title("$\Delta$Distance Along Arc")
            plt.xlabel("x (px)")
            plt.ylabel("$\Delta$Distance along arc (px)")
            plt.plot(points[:-1,1],np.diff(distance) * self.config["convert"]["mmtopx"],'g.')
            plt.savefig("%sInstrument-%04d-Delta-Distances%s" % (self.config["Dirs"]["Partials"],lenslet_num,self.config["plot_format"]))
            plt.clf()
        
        xs,ys = points.T
        return xs,ys,wl
    
    
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
        xs += (self.config["image_size"]["mm"]/2)
        ys += (self.config["image_size"]["mm"]/2)
        
        lams *= 1e-6 # Convert wavelength to SI units (m)
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
        
        # Progress bar for lenslet creation and validation
        PBar = arpytools.progressbar.ProgressBar(color="green")
        finished = 0.0
        total = len(self.lenslets)
        self.log.useConsole(False)
        PBar.render(0,"L:%4d %4d/%-4d" % (0,finished,total))
        
        # Variables for lenslet use
        self.lensletObjects = {}
        FileName = self.config["Dirs"]["Partials"] + "Lenslets-raw" + ".dat"
        with open(FileName,'w') as stream:
            for idx in self.lenslets:
                select = idx == ix
                aLenslet = Lenslet(xs[select],ys[select],xpix[select],ypix[select],p1[select],p2[select],lams[select],idx,self.config)
                if aLenslet.valid():
                    self.lensletObjects[idx] = aLenslet
                stream.write(aLenslet.introspect())
                progress = int((finished/float(total)) * 100)
                finished += 1
                PBar.render(progress,"L:%4d %4d/%-4d" % (idx,finished,total))
        PBar.render(100,"L:%4d %4d/%-4d" % (idx,total,total))
        self.lenslets = self.lensletObjects.keys()
        self.log.useConsole(True)
        
    def plot_lenslet_data(self):
        """Outputs the lenslet data"""
        self.log.info("Generating Lenslet Plots")
        plt.clf()
        FileName = "%(dir)sLenslet-xy%(fmt)s" % { 'dir' : self.config["Dirs"]["Partials"], 'fmt':self.config["plot_format"]}
        for lenslet in self.lensletObjects.values():
            plt.plot(lenslet.xs,lenslet.ys,linestyle='-')
        plt.title("Lenslet x-y positions")
        plt.savefig(FileName)
        plt.clf()
        FileName = "%(dir)sLenslet-pxy%(fmt)s" % { 'dir' : self.config["Dirs"]["Partials"], 'fmt':self.config["plot_format"]}
        for lenslet in self.lensletObjects.values():
            x,y = lenslet.ps.T
            plt.plot(x,y,marker='.')
        plt.title("Lenslet p-xy positions")
        plt.savefig(FileName)
        plt.clf()
    
    def loadOpticsData(self,laspec,dispspec):
        """Loads an optical conversion based on the lenslet array spec and dispersion spec files provided. This wrapper function handles both file loading functions. See `loadDispersionData` and `loadLensletData`."""
        self.log.info("Loading Optics Data from Spec Files")
        self.loadDispersionData(dispspec)
        self.loadLensletData(laspec)
    
    
    def setupLenslets(self):
        """Establish the list of lenslets for use in the system"""        
        if "start" in self.config["Lenslets"]:
            self.lenslets = self.lenslets[self.config["Lenslets"]["start"]:]
        if "number" in self.config["Lenslets"]:
            self.lenslets = self.lenslets[:self.config["Lenslets"]["number"]]
        self.total = len(self.lenslets)
            
    # Image Tracing Function
    def cache_images(self,get_spec):
        """Retrieve dense images for everything"""
        PBar = arpytools.progressbar.ProgressBar(color="blue")
        finished = 0.0
        total = len(self.lenslets)
        
        self.log.info("Generating Subimages for %d lenslets" % total)
        
        self.log.toggleConsole(False)
        PBar.render(0,"L:%4d %4d/%4d" % (0,finished,total))
        
        for index in self.lenslets:
            self.cache_sed_subimage(index,get_spec(index),write=self.config["Cache"])
            finished += 1.0
            progress = int((finished/float(total)) * 100)
            PBar.render(progress,"L:%4d %4d/%-4d" % (index,finished,total))
        
        self.log.toggleConsole(True)
        
    
    def get_dense_image(self,lenslet_num,spectrum):
        """docstring for get_dense_image"""
        x,y,flux,wl,xsize,ysize,corner = self.get_trace(lenslet_num,spectrum)
        
        img = np.zeros((xsize,ysize))
        
        img[x,y] = flux
        self.log.debug("Placing flux (shape: %s ) for spectrum %4d into a sub-image (shape: %s)."
            % (flux.shape,lenslet_num,img.shape))
        
        return img,corner
        
    
    def get_trace(self,lenslet_num,spectrum):
        """This function returns a dense image array, with flux placed into single pixels. The image returned is too dense for use on the final CCD, and must be also be convolved with the telescope image and spectral PSF. The image, and the numpy array first corner is returned by this function.
        
        .. Note:: The source code of this function is heavily documented."""
        lenslet = self.lensletObjects[lenslet_num]
        dx,dy,wl = lenslet.dxs,lenslet.dys,lenslet.dwl
        points = np.array([dx,dy]).T
        deltawl = np.diff(wl)
        # Some notes about this function:
        #  points:  These points are in normal CCD Pixel units, but are not integers,
        #           as they correspond to dense pixels on the overdense images
        #  wl:      A list of wavelengths that correspond to the points above
        
        # Create the spectra, by calling the spectra with wavelength
        # (in meters, note the conversion from wl, which was originally in microns)
        # Then use the deltawl to get the true amount of flux in each pixel
        self.log.debug(npArrayInfo(wl,"Spectrum wl Request"))
        
        WLS = wl*1e-6
        RS = np.diff(WLS) 
        WLS = WLS[:-1]
        RS = WLS/RS
        
        
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
            with open(self.config["Dirs"]["Partials"]+"Instrument-Offsets.dat",'a') as handle:
                np.savetxt(handle,offset)
        corner -= np.array([-self.config["padding"],self.config["padding"]])
        
        x += offset[0]
        y += offset[1]
        
        # Create our sub-image, using the x and y width of the spectrum, plus 2 padding widths.
        # Padding is specified in full-size pixels to ensure that the final image is an integer
        # number of full-size pixels across.
        xsize = xdist+2*self.config["padding"]*self.config["density"]
        ysize = ydist+2*self.config["padding"]*self.config["density"]
        
        self.log.debug("Scailing by %g" % self.config["gain"])
        radiance = spectrum(wavelengths=WLS,resolution=RS) 
        radiance *= self.config["gain"]
        self.log.debug(npArrayInfo(radiance,"Generated Spectrum"))
        self.log.debug(npArrayInfo(deltawl,"DeltaWL Rescaling"))
        flux = radiance[1,:] * deltawl *1e-6
        self.log.debug(npArrayInfo(flux,"Final Flux"))
        
        # This debugging area generates plots for us.
        if self.log.getEffectiveLevel() <= logging.DEBUG and self.config["plot"]:
            self.log.debug("Generating Plots for Spectra...")
            plt.clf()
            plt.plot(WLS,radiance[1,:],"b.")
            plt.title("Retrieved Spectra")
            plt.xlabel("Wavelength ($\mu m$)")
            plt.ylabel("Radiance (Units undefined)")
            plt.savefig("%sInstrument-%04d-Radiance%s" % (self.config["Dirs"]["Partials"],lenslet,self.config["plot_format"]))
            plt.clf()
            plt.clf()
            plt.plot(wl[:-1]*1e-6,flux,"b.")
            plt.title("Generated, Fluxed Spectra")
            plt.xlabel("Wavelength ($\mu m$)")
            plt.ylabel("Flux (Units undefined)")
            plt.savefig("%sInstrument-%04d-Flux%s" % (self.config["Dirs"]["Partials"],lenslet,self.config["plot_format"]))
            plt.clf()
            plt.plot(wl[:-1]*1e-6,deltawl*1e-6,"g.")
            plt.title("$\Delta\lambda$ for each pixel")
            plt.xlabel("Wavelength ($\mu m$)")
            plt.ylabel("$\Delta\lambda$ per pixel")
            plt.savefig("%sInstrument-%04d-DeltaWL%s" % (self.config["Dirs"]["Partials"],lenslet,self.config["plot_format"]))
            plt.clf()
            plt.semilogy(WLS,RS,"g.")
            plt.title("$R = \\frac{\lambda}{\Delta\lambda}$ for each pixel")
            plt.xlabel("Wavelength ($\mu m$)")
            plt.ylabel("Resolution $R = \\frac{\Delta\lambda}{\lambda}$ per pixel")
            plt.savefig("%sInstrument-%04d-Resolution%s" % (self.config["Dirs"]["Partials"],lenslet,self.config["plot_format"]))
            plt.clf()
        
        
        
        img = np.zeros((xsize,ysize))
        
        # imshape = self.Caches.get("CONV").shape
        #        
        #        for x,y,flux in zip(x,y,flux,WLS):
        #            tinyImage = np.zeros(imshape)
        #            corner = [ x + imshape[0]/2.0, y + imshape[0]/2.0]
        #            self.log.debug("Corner of super-small image: %s" % corner)
        #            tinyImage[x,y] = flux
        #            PSF = self.get_psf(wl)
        #            tinyImage = sp.convolve(tinyImage,PSF,mode="same")
        
        
        # Place the spectrum into the sub-image
 
        
        if self.log.getEffectiveLevel() <= logging.DEBUG:
            np.savetxt(self.config["Dirs"]["Partials"]+"Instrument-Subimage-Values.dat",np.array([x,y,wl[:-1],deltawl,flux]).T)
        
        return x,y,flux,wl[:-1],xsize,ysize,corner
        
    
    
    # Image Convolution
    def get_sub_image(self,lenslet,spectrum,fast=False):
        """Returns a sub-image for a given lenslet, which has been binned and convolved with the appropriate telescope and PSF images. There are two modes for this fucntion: fast, which uses a single convolution, and slow, which does each convolution separately and returns intermediate stage images for diagnostics."""
        # This function gets the dense sub-image with the spectrum placed in single dense pixels
        img,corner = self.get_dense_image(lenslet,spectrum)
        self.log.debug("Retrieved Dense Image for %4d" % lenslet)
        
        if fast:
            # Convolve with the PSF and Telescope Image simultaneously
            img2 = sp.signal.convolve(img,self.Caches.get("CONV"),mode='same')
            self.log.debug("Convolved Dense Image with PSF and Telescope for %4d" % lenslet)
            # Bin the image back down to the final pixel size
            small = bin(img2,self.config["density"]).astype(np.int16)
            self.log.debug("Binned Dense Image to Actual Size for %4d" % lenslet)
            return small,corner
        else:
            # Convolve this spectrum with an appropriate image of the telescope
            img_tel = sp.signal.convolve(img,self.Caches.get("TEL"),mode='same')
            self.log.debug("Convolved Dense Image with Telescope for %4d" % lenslet)
            # Convolve again with an appropriate PSF
            img_tel_psf = sp.signal.convolve(img_tel,self.Caches.get("PSF"),mode='same')
            self.log.debug("Convolved Dense Image with PSF for %4d" % lenslet)
            # Bin the image back down to the final pixel size
            small = bin(img_tel_psf,self.config["density"]).astype(np.int16)
            self.log.debug("Binned Dense Image for %4d" % lenslet)
            return small,corner,(img,img_tel,img_tel_psf)
    
    
    # Image Caching
    def cache_sed_subimage(self,lenslet,spectrum,write=False,do_return=False):
        """Generates a sub image, and saves that result to this object. Should be thread-safe."""
        if self.debug:
            small, corner, steps = self.get_sub_image(lenslet,spectrum)
        else:
            small, corner = self.get_sub_image(lenslet,spectrum,fast=True)
        
        self.log.debug("Retrieved SED Subimage for lenslet %d" % lenslet)
        
        label = "SUBIMG%d" % lenslet
        
        self.save(small,label)
        frame = self.frame()
        frame.lensletNumber = lenslet
        frame.corner = corner
        frame.spectrum = spectrum.label
        frame.config = hash(str(self.config))
        

        if self.debug and self.config["plot"]:
            Stages = ["Raw Image","Convolved with Telescope","Convolved with Telescope and PSF"]
            StagesF = ["Raw","Tel","PSF"]
            for i,step in enumerate(steps):
                self.save(step,"%04d-Intermediate-%d: %s" % (lenslet,i,Stages[i]))
                frame = self.frame()
                frame.lensletNumber = lenslet
                frame.corner = corner
                frame.spectrum = spectrum.label
                frame.config = hash(str(self.config))
                frame.header.update(dict(SEDstg=Stages[i]))
                plt.imshow(step)
                plt.title("%s Image Generation Steps for Lenslet %4d" % (Stages[i],lenslet))
                plt.savefig("%sInstrument-%04d-Subimage-%s%s" % (self.config["Dirs"]["Partials"],lenslet,StagesF[i],self.config["plot_format"]))
                plt.clf()
        
        # We only write the sub-image if the function is called to write sub images
        if write:
            toSave = [label]
            if self.debug and self.config["plot"]:
                toSave += ["%04d-Intermediate-%d: %s" % (lenslet,i,Stages[i]) for i in range(len(steps))]
            self.write("%sSubimage-%4d%s" % (self.config["Dirs"]["Caches"],lenslet,".fits"),toSave,label,clobber=True)
            self.remove(label)
            if self.debug and self.config["plot"]:
                for i,step in enumerate(steps):
                    self.remove("%04d-Intermediate-%d: %s" % (lenslet,i,Stages[i]))
                    
        
        if do_return:
            return small, corner
        return
    
    
    # Placement Functions
    def place_cached_sed(self,lenslet,label,dlabel,fromfile=False):
        """Places a cached SED Subimage"""
        slabel = "SUBIMG%d" % lenslet
        subframe = False
        if fromfile:
            try:
                filename = "%sSubimage-%4d%s" % (self.config["Dirs"]["Caches"],lenslet,".fits")
                slabel = self.read(filename)[0]
                subframe = self.frame(slabel)
            except IOError as e:
                self.log.debug("Could not load spectrum %s from file" % filename)
                raise SEDLimits(str(e))
        else:
            try:
                subframe = self.frame(slabel)
            except KeyError as e:
                raise SEDLimits(str(e))
        subimg = subframe()
        mlenslet = subframe.lensletNumber
        mcorner = subframe.corner
        if mlenslet != lenslet:
            raise ValueError("Lenslet Number Mismatch %d:%d for state %s in %s" % (lenslet,mlenslet,slabel,self))
        self.remove(slabel)
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
        # self.log.debug("Placed a new image %s into %s at corner %s" % (label,dlabel,corner))
        self.remove(label)
        
    
    
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
    



