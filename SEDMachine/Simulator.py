#!/usr/bin/env python
# 
#  Simulator.py
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2012-02-08.
#  Copyright 2012 Alexander Rudy. All rights reserved.
#  Version 0.3
# 

import numpy as np
# import pyfits as pf
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.signal
import scipy.interpolate
import yaml

import shapely as sh
import shapely.geometry

import arpytools.progressbar

import os
import logging,logging.handlers
import time
import copy
import collections

import AstroObject
from AstroObject.AstroSimulator import Simulator
from AstroObject.AstroCache import *
from AstroObject.AstroSpectra import SpectraObject
from AstroObject.AstroImage import ImageObject,ImageFrame
from AstroObject.AnalyticSpectra import BlackBodySpectrum, AnalyticSpectrum, FlatSpectrum, ResampledSpectrum
from AstroObject.Utilities import *

from Lenslet import *

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


class SEDSimulator(Simulator,ImageObject):
    """A simulator for the SED Machine"""
    def __init__(self):
        super(SEDSimulator, self).__init__(name="SEDMachine")
        self.debug = False
        self.dataClasses = [SubImage]
        self.lenslets = []
        self.config.merge(self.basics)
        self.config.merge({"Instrument":self.instrument,"Caches":self.caches,"Source":self.source})
        self.config.setFile("Main")
        self.setup_stages()

    
    basics = {
        "Cache": True,
        "Debug": False,
        "Plot": False,
        'plot_format': '.pdf', 
        "Output": {
            "Label": "Generated",
            "Format": "fits",
        },
        "Configurations": {
            "Instrument" : "SED.instrument.config.yaml",
            "Source" : "SED.source.config.yaml",
            "Main" : "SED.main.config.yaml",
        },
        "Dirs": {
            "Logs": "Logs",
            "Partials" : "Partials",
            "Caches" : "Caches",
            "Images" : "Images",
        },
        "Lenslets" : {},
        'logging': {
            'console': {
                'level': False, 
                'enable': True, 
                'format': '... %(levelname)-8s %(message)s'
            },
            'growl' : {
                'name' : "SED Machine Simulator",
            },
            'file': {
                'format': '%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s',
                'enable': True, 
                'filename': 'SEDMachine'
            },
        },
    }
    
    caches = {
        'Telescope' : "SED.Instrument.tel.npy",
        'PSF' : "SED.Instrument.psf.npy",
        'CONV' : "SED.Instrument.conv.npy",
        'config' : "SED.Instrument.config.yaml",
    }
       
    instrument = {
        'files': {
            'dispersion': 'Data/dispersion_12-10-2011.txt',
            'encircledenergy': 'Data/encircled_energy_4nov11.TXT',
            'lenslets': 'Data/xy_17nov2011_v57.TXT',
        },
        'dark': 20, 
        'convert': {
            'pxtomm': 0.0135 }, 
        'density': 5, 
        'tel_obsc': {'px': 0.2}, 
        'plot': False, 
        'ccd_size': {'px': 2048}, 
        'padding': 5, 
        'psf_stdev': {'px': 1.0}, 
        'bias': 20, 
        'psf_size': { 'px': 0}, 
        'image_size': { 'mm': 40.0}, 
        'tel_radii': { 'px': 1.2}, 
        'exposure': 120, 
        'gain': 10000000000.0,
    }
    
    source = {
        'Filename' : "Data/SNIa.R1000.dat",
        'PreAmp' : 1.0,
    }
    
    def setup_stages(self):
        """Sets up all simulator stages"""
        self.registerStage(self.setup_caches,"setup_caches",help=False,description="Setting up caches")
        self.registerStage(self.setup_configuration,"setup_config",help=False,description="Setting up dynamic configuration")
        self.registerStage(self.setup_lenslets,"setup_lenslets",help=False,description="Setting up lenslets",dependencies=["setup_config"])
        self.registerStage(self.setup_blank,"setup_blank",help=False,description="Creating blank image")
        self.registerStage(self.setup_source,"setup_source",help=False,description="Creating Spectrum Objects")
        self.registerStage(None,"setup",dependencies=["setup_caches","setup_lenslets","setup_blank","setup_source"],help="System Setup",description="Set up simulator")
        self.registerStage(self.plot_lenslet_data,"plot_lenslets",help="Plot Lenslets",description="Plotting Lenslet Positions",include=False,dependencies=["setup_lenslets"])
        self.registerStage(self.lenslet_dispersion,"dispersion",help="Calculate dispersion",description="Calculating Dispersion for each Lenslet",dependencies=["setup_lenslets","setup_caches"])
        self.registerStage(self.lenslet_trace,"trace",help="Trace Lenslets",description="Tracing Lenslet Dispersion",dependencies=["dispersion","setup_caches","setup_source"])
        self.registerStage(self.lenslet_place,"place",help="Place Subimages",description="Placing Lenslet Spectra",dependencies=["trace"])
        self.registerStage(self.lenslet_write,"write",help="Write Subimages",description="Writing Subimages to Cache Files",dependencies=["place"])
    
    def setup_caches(self):
        """Register all of the cache objects and types"""
        self.Caches.registerNPY("TEL",self.get_tel_kern,filename=self.config["Caches"]["Telescope"])
        self.Caches.registerNPY("PSF",self.get_psf_kern,filename=self.config["Caches"]["PSF"])
        self.Caches.registerNPY("CONV",lambda : sp.signal.convolve(self.Caches.get("PSF"),self.Caches.get("TEL"),mode='same'),filename=self.config["Caches"]["CONV"])
        self.Caches.registerCustom("CONFIG",kind=YAMLCache,generate=self.config.extract,filename=self.config["Caches"]["config"])
        
        if "clear_cache" in self.options and self.options["clear_cache"]:
            self.Caches.clear()
        if "cache" in self.options and not self.options["clear_cache"]:
            self.Caches.disable()
        
        if self.Caches.check(master="CONFIG"):
            self.log.info("Caches appear to be out of date, regenerating")
        self.Caches.get("CONFIG")
        
    def setup_lenslets(self):
       """This function loads data about lenslet positions, and thier dispersion through the prism. The data are original produced by Zeemax. This function reads the Zeemax data directly and then cleans the data in certain ways, preparing it for use later in the system.
        
       ..Note:: The source of this function is well documented.
        
       ..Note:: This function does not store variables neatly. As such, it has no built-in caching system.
       """
       # Load Lenslet Specification File
       ix, p1, p2, lams, xs, ys = np.genfromtxt(self.config["Instrument"]["files"]["lenslets"],skip_header=1).T
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
        
       lams *= 1e-6 # Convert wavelength to SI units (m)
       # Find the xs and ys that are not within 0.1 mm of the edge of the detector...
       ok = (xs > 0.1) & (xs < self.config["Instrument"]["image_size"]["mm"]-0.1) & (ys > 0.1) & (ys < self.config["Instrument"]["image_size"]["mm"]-0.1)
       ix, p1, p2, lams, xs, ys = ix[ok], p1[ok], p2[ok], lams[ok], xs[ok], ys[ok]
       # We remove these positions because we don't want to generate spectra for them.
        
       # This simply generates a list of all of the lenslets
       self.lensletIndex = np.unique(ix)
        
       # Convert the xs and ys to pixel positions
       xpix = np.round(xs * self.config["Instrument"]["convert"]["mmtopx"],0).astype(np.int)
       ypix = np.round(ys * self.config["Instrument"]["convert"]["mmtopx"],0).astype(np.int)
        
       # Determine the center of the whole system by finding the x position that is closest to 0,0 in pupil position
       cntix = np.argmin(p1**2 + p2**2)
       self.center = (xs[cntix] * self.config["Instrument"]["convert"]["mmtopx"], ys[cntix] * self.config["Instrument"]["convert"]["mmtopx"])
        
       # Progress bar for lenslet creation and validation
       PBar = arpytools.progressbar.ProgressBar(color="green")
       finished = 0.0
       total = len(self.lensletIndex)
       self.log.useConsole(False)
       PBar.render(0,"L:%4s %4d/%-4d" % ("",finished,total))
        
       # Variables for lenslet use
       self.lenslets = {}
       FileName = self.config["Dirs"]["Partials"] + "Lenslets-raw" + ".dat"
       with open(FileName,'w') as stream:
           for idx in self.lensletIndex:
               select = idx == ix
               lenslet = Lenslet(xs[select],ys[select],xpix[select],ypix[select],p1[select],p2[select],lams[select],idx,self.config,self.Caches)
               if lenslet.valid():
                   self.lenslets[idx] = lenslet
                   stream.write(lenslet.introspect())
               progress = int((finished/float(total)) * 100)
               finished += 1
               PBar.render(progress,"L:%4d %4d/%-4d" % (idx,finished,total))
       PBar.render(100,"L:%4s %4d/%-4d" % ("Done",total,total))
       self.lensletIndex = self.lenslets.keys()
       self.log.useConsole(True)
       if "start" in self.config["Lenslets"]:
           self.lensletIndex = self.lensletIndex[self.config["Lenslets"]["start"]:]
       if "number" in self.config["Lenslets"]:
           self.lensletIndex = self.lensletIndex[:self.config["Lenslets"]["number"]]
       self.total = len(self.lensletIndex)
    
    def setup_blank(self):
        """Establish a blank Image"""
        self.save(np.zeros((self.config["Instrument"]["image_size"]["px"],self.config["Instrument"]["image_size"]["px"])).astype(np.int16),"Blank")

    def setup_source(self):
        """Sets up a uniform source file based spectrum"""
        WL,Flux = np.genfromtxt(self.config["Source"]["Filename"]).T
        WL *= 1e-10
        Flux /= np.max(Flux) 
        self.R_Spectrum = ResampledSpectrum(np.array([WL,Flux]),self.config["Source"]["Filename"])
        self.Spectrum = self.R_Spectrum * self.config["Source"]["PreAmp"]
    

    def lenslet_dispersion(self):
        """Calculate the dispersion for each lenslet"""
        self.map_over_lenslets(lambda l: l.find_dispersion(),exceptions=None,color="blue")
        
    def lenslet_trace(self):
        """Trace out each lenslet"""
        self.map_over_lenslets(lambda l: l.get_trace(self.Spectrum),exceptions=None,color="blue")

    def lenslet_place(self):
        """Place each spectrum into the subimage"""
        self.map_over_lenslets(lambda l: l.place_trace(self.get_psf),color="yellow")
        
    def lenslet_write(self):
        """Write each lenslet into the final image"""
        self.map_over_lenslets(lambda l: l.write("%sSubimage-%4d%s" % (self.config["Dirs"]["Caches"],lenslet,".fits")),color="red")
        
    def lenslet_delete(self):
        """Clear out lenslet subimages for memory purposes"""
        self.map_over_lenslets(lambda l: del l,color="red")
    
    #######################
    ## DEBUGGING METHODS ##
    #######################
    
    def plot_lenslet_data(self):
        """Outputs the lenslet data"""
        plt.figure()
        plt.clf()
        
        FileName = "%(dir)sLenslet-xy%(fmt)s" % { 'dir' : self.config["Dirs"]["Partials"], 'fmt':self.config["plot_format"]}
        for lenslet in self.lenslets.values():
            plt.plot(lenslet.xs,lenslet.ys,linestyle='-')
        plt.title("Lenslet x-y positions")
        plt.savefig(FileName)
        
        plt.clf()
        
        FileName = "%(dir)sLenslet-pxy%(fmt)s" % { 'dir' : self.config["Dirs"]["Partials"], 'fmt':self.config["plot_format"]}
        for lenslet in self.lenslets.values():
            x,y = lenslet.ps.T
            plt.plot(x,y,marker='.')
        plt.title("Lenslet p-xy positions")
        plt.savefig(FileName)
        
        plt.clf()
    
    def map_over_lenslets(self,function,exceptions=True,color="green"):
        """Maps a given function to operate on each lenslet, and displays a progress bar along the way."""
        if exceptions == True:
            exceptions = Exception
        
        self.progress = 0.0
        self.errors = 0
        self.bar = arpytools.progressbar.ProgressBar(color=color)
        self.log.useConsole(False)
        self.bar.render(0,"L:%4s %4d/%-4d" % ("",self.progress,self.total))
        for lenslet in self.lenslets.values():
            self.bar.render(int(self.progress/self.total * 100),"L:%4d %4d/%-4d" % (lenslet.num,self.progress,self.total))
            try:
                function(lenslet)
            except exceptions as e:
                self.log.error("Caught %s in Lenslet %d" % (e.__class__.__name__,lenslet.num))
                self.log.error(str(e))
                self.errors += 1
            finally:
                self.progress += 1.0
                self.bar.render(int(self.progress/self.total * 100),"L:%4d %4d/%-4d" % (lenslet.num,self.progress,self.total))
        self.bar.render(100,"L:%4s %4d/%-4d" % ("Done",self.progress,self.total))
        self.log.useConsole(True)
        if self.progress != self.total:
            self.log.warning("Progress and Total are different at end of loop: %d != %d" % (self.progress,self.total))
        if self.errors != 0:
            self.log.warning("Trapped %d errors" % self.errors)
        
        
    ###################
    ## Image KERNELS ##
    ###################
    
    def get_psf(self,wavelength):
        """Return a PSF for a given wavelength in the system"""
        return self.Caches.get("PSF")
    
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
        TELIMG = self.circle_kern( self.config["Instrument"]["tel_radii"]["px"] * self.config["density"] )
        center = self.circle_kern( self.config["Instrument"]["tel_obsc"]["px"] * self.config["density"] ,
            self.config["Instrument"]["tel_radii"]["px"] * self.config["density"] , False )
        TELIMG -= center
        TELIMG = TELIMG / np.sum(TELIMG)
        self.log.debug("Generated a Telescpe Kernel with shape %s" % (str(TELIMG.shape)))
        return TELIMG
    
    def get_psf_kern(self):
        """Returns the PSF Kernel. The function first tries to read the encircled energy file. In this case, if a `psf_size` is set in the instrument configuration, this value will be used to truncate the size of the encircled energy function. If the encircled energy function cannot be loaded, the system will fall back on to a gaussian psf as configured by the instrument."""
        if self.config["Instrument"]["psf_size"]["px"] != 0:
            size = self.config["Instrument"]["psf_size"]["px"] * self.config["density"]
            truncate = True
        else:
            size = 0
            truncate = False
        try:
            PSFIMG = self.psf_kern( self.config["Instrument"]["files"]["encircledenergy"],size,truncate)
        except IOError as e:
            self.log.warning("Could not access encircled energy file: %s" % e)
            PSFIMG = self.gauss_kern( (self.config["Instrument"]["psf_stdev"]["px"] * self.config["density"]) )
        else:
            self.log.debug("Loaded Encircled Energy from %s" % self.config["Instrument"]["files"]["encircledenergy"])
        return PSFIMG
        
    def setup_configuration(self):
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
    
SIM = SEDSimulator()
SIM.run()  