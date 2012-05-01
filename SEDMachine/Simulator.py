#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  Simulator.py
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2012-02-08.
#  Copyright 2012 Alexander Rudy. All rights reserved.
#  Version 0.3.9
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

import os
import logging,logging.handlers
import time
import copy
import collections
import gc

from pkg_resources import resource_filename

import AstroObject
from AstroObject.AstroSimulator import *
from AstroObject.AstroCache import *
from AstroObject.AstroConfig import *
from AstroObject.AstroSpectra import SpectraStack,SpectraFrame
from AstroObject.AstroImage import ImageStack,ImageFrame
from AstroObject.AnalyticSpectra import BlackBodySpectrum,GaussianSpectrum, AnalyticSpectrum, FlatSpectrum, InterpolatedSpectrum, UnitarySpectrum
from AstroObject.Utilities import *
from AstroObject.mpl_utilities import *

from version import version as versionstr
from Objects import *


class SEDSimulator(Simulator,ImageStack):
    """A simulator implementation for the SED Machine.
    
    This simulator is based on :class:`AstroObject.AstroSimulator.Simulator`. It is designed to run as a series of dependent stages. The simulator is first setup, with basic data structures created in the constructor, and then simulator stages registered in :meth:`setup_stages`."""
    def __init__(self):
        super(SEDSimulator, self).__init__(name="SEDMachine",version=versionstr)
        self.debug = False
        self.mapping = False
        self.dataClasses = [SubImage]
        self.lenslets = {}
        self.ellipses = {}
        self.qe = SpectraStack(dataClasses=[AnalyticSpectrum,SpectraFrame])
        self.qe.save(FlatSpectrum(0.0))
        self.spectra =  SpectraStack(dataClasses=[AnalyticSpectrum,SpectraFrame])
        self.spectra.save(FlatSpectrum(0.0))
        self.sky =  SpectraStack(dataClasses=[AnalyticSpectrum,SpectraFrame])
        self.sky.save(FlatSpectrum(0.0))
        self.astrologger = logging.getLogger("AstroObject")
        self.config.load(resource_filename(__name__,"SED.main.config.default.yaml"))
        self.config.setFile("Main")
        self.setup_stages()
    
    def setup_stages(self):
        """Registers all of the availalbe simulator stages. For basic command help on the options registered here, use::
            
            $ SEDMsim --help
            
        To get a list of all available stages, use::
            
            $ SEDMsim --stages
        
        This function is **ONLY** used to register stages and configuration options with the command-line parser. It is called by :class:`SEDSimulator` during construction.
        """
        self.registerConfigOpts("D",{"Lenslets":{"start":2100,"number":50},"Debug":True,"Output":{"Label":"DFlag",},},help="Debug, Limit lenslets (50,start from 2100)")
        self.registerConfigOpts("S",{"Lenslets":{"start":2100,"number":5},"Debug":True,"Output":{"Label":"SFlag",},},help="Debug, Limit lenslets (5,start from 2100)")
        self.registerConfigOpts("T",{"Lenslets":{"start":2100,"number":50},"Debug":False,"Output":{"Label":"TFlag",},},help="Limit lenslets (50,start from 2100)")
        self.registerConfigOpts("M",{"Lenslets":{"start":1000,"number":500},"Debug":True,"Output":{"Label":"MFlag",},},help="Limit lenslets (500,start from 1000)")
        self.registerConfigOpts("N",{"Lenslets":{"start":1000,"number":500},"Debug":False,"Output":{"Label":"NFlag",},},help="Limit lenslets (500,start from 1000)")
        self.registerConfigOpts("A",{"Lenslets":{"start":2100,"number":1},"Debug":True,"Output":{"Label":"AFlag",},},help="Debug, Single lenslets (start from 2100)")
        self.registerConfigOpts("C",{"Lenslets":{"position":{"x":0.0,"y":0.0},"radius":0.01},"Debug":True,"Output":{"Label":"CFlag",},},help="Debug, Central Lenslets")
        
        # SETUP Stages
        self.registerStage(self.setup_caches,"setup-caches")
        self.registerStage(self.setup_configuration,"setup-config")
        self.registerStage(self.setup_constants,"setup-constants")
        self.registerStage(self.setup_cameras,"setup-cameras")
        self.registerStage(self.setup_lenslets,"setup-lenslets")
        self.registerStage(self.setup_hexagons,"setup-hexagons")
        self.registerStage(self.setup_blank,"setup-blank")
        self.registerStage(self.setup_dummy_blank,"setup-blank-d")

        self.registerStage(self.setup_simple_source,"setup-source-simple")
        self.registerStage(None,"simple-source",help="Use a simple, centered source object",description="Replacing default source with a simple one",dependencies=["setup-source-simple"])
        self.registerStage(None,"setup-source",help=False,description="Creating source spectrum objects",dependencies=["setup-config","setup-constants","simple-source"])
        self.registerStage(self.setup_source_pixels,"setup-source-pixels")
        self.registerStage(self.setup_noise,"setup-noise")
        self.registerStage(self.setup_sky,"setup-sky")
        self.registerStage(self.setup_line_list,"setup-lines")
        self.registerStage(self.setup_scatter,"setup-scatter")
        # Setup Macro
        self.registerStage(None,"setup",help="System Setup",description="Set up simulator",
            dependencies=["setup-hexagons","setup-blank","setup-source","setup-noise","setup-constants","setup-sky","setup-cameras","setup-lines","setup-scatter"],
            )
        
        
        self.registerStage(self.geometric_resample,"geometric-resample")
        
        # Apply spectral properties
        self.registerStage(self.apply_sky,"apply-sky")
        self.registerStage(self.apply_qe,"apply-qe")
        self.registerStage(self.apply_atmosphere,"apply-atmosphere")

        # Adjust spectra
        self.registerStage(self.sky_source,"sky-source")
        self.registerStage(self.flat_source,"flat-source")
        self.registerStage(self.line_source,"line-source")
        
        
        # Dispersion functions
        self.registerStage(self.lenslet_dispersion,"dispersion")
        self.registerStage(self.lenslet_trace,"trace")
        self.registerStage(self.lenslet_place,"place")
        
        # Merge images back together
        self.registerStage(self.image_merge,"merge-cached")
        self.registerStage(None,"merge",help="Merge Subimages",description="Merging subimage into master image",dependencies=["place","merge-cached"])
        
        # Final Image work
        self.registerStage(self.ccd_crop,"crop")
        self.registerStage(self.apply_scatter,"add-scatter")
        self.registerStage(self.apply_noise,"add-noise")
        self.registerStage(self.transpose,"transpose")
        self.registerStage(self.save_file,"save")
        
        # Alternative work macros
        self.registerStage(None,"cached-only",help="Use cached subimages to construct final image",description="Building image from caches",dependencies=["merge-cached","crop","add-noise","add-scatter","transpose","save"])
        
        
        # Plotting geometry functions
        self.registerStage(self.plot_kernel_partials,"plot-kernel")
        self.registerStage(self.plot_hexagons,"plot-hexagons")
        self.registerStage(self.plot_invalid_hexagons,"plot-bad-hexagons")
        self.registerStage(self.plot_pixels,"plot-pixels")
        self.registerStage(self.plot_invalid_pixels,"plot-bad-pixels")
        self.registerStage(self.plot_geometry,"plot-p-geometry")
        self.registerStage(self.write_resample,"write-resample")
        self.registerStage(self.plot_resample,"plot-resample")
        self.registerStage(None,"plot-geometry",help="Do all geometry plots",description="Plotting geometries",dependencies=["plot-hexagons","plot-invalid-hexagons","plot-pixels","plot-bad-pixels","plot-p-geometry","plot-resample"])
        
        # Spectrum Plotting Functions
        self.registerStage(self.compare_methods,"plot-spectrum-tests")
        self.registerStage(self.plot_original_calibration,"plot-cal-o")
        self.registerStage(self.plot_original_source,"plot-source-o")
        self.registerStage(self.plot_source,"plot-source")
        self.registerStage(self.plot_sky_original,"plot-sky-o")
        self.registerStage(self.plot_sky,"plot-sky")
        self.registerStage(self.plot_qe,"plot-qe")
        
        self.registerStage(None,"plot-spectra",help=False,dependencies=["plot-qe","plot-sky","plot-source","plot-cal-o"])
        
        
        # Dispersion plotting functions
        self.registerStage(self.plot_ellipses,"plot-lenslet-es")
        self.registerStage(self.plot_rotation,"plot-lenslet-rs")
        
        
        self.registerStage(self.plot_lenslet_data,"plot-lenslet-xy")
        self.registerStage(self.plot_dispersion_data,"plot-dispersion")
        self.registerStage(self.plot_trace_data,"plot-trace")
        self.registerStage(self.plot_spectrum_data,"plot-spectrum")
        self.registerStage(None,"plot-lenslets",help=False,description="Plotting data about each lenslet",dependencies=["plot-dispersion","plot-trace","plot-spectrum"])
        
        
        self.registerStage(None,"plot",help="Create all plots",description="Plotting everything",dependencies=["plot-lenslet-xy","plot-lenslets","plot-sky","plot-qe","plot-source","plot-p-geometry","plot-hexagons","plot-pixels","plot-spectrum-tests","plot-kernel"])
        
    
    @description("Setting up Caches")
    def setup_caches(self):
        """Establish Cache objects for image kernels. The Caching algorithm doesn't quite work properly yet, so this method is really superfluous.
        
        **Command Name:** ``*setup-caches``
        """
        for key in self.config["Caches"]:
            self.config["Caches"][key] = "%s/%s" % (self.config["Dirs.Caches"],self.config["Caches"][key])
        self.Caches["TEL"] = NumpyCache(self.get_tel_kern,filename=self.config["Caches.Telescope"])
        self.Caches["PSF"] = NumpyCache(self.get_psf_kern,filename=self.config["Caches.PSF"])
        self.Caches["CONV"] = NumpyCache(lambda : sp.signal.convolve(self.Caches["PSF"],self.Caches["TEL"],mode='same'),filename=self.config["Caches.CONV"])
        
        if "clear_cache" in self.config["Options"] and self.config["Options"]["clear_cache"]:
            self.Caches.flag('enabled',False)
        if "cache" in self.config["Options"] and not self.config["Options"]["cache"]:
            self.Caches.flag('saving',False)
    
    
    
    @description("Setting up operational constants")
    @depends("setup-caches")
    def setup_constants(self):
        """Generate the physical constants required for simple calculations.
        
        **Command Name:** ``*setup-constants``
        
        
        **Constants Generated**
        
        :var hc: Plank's constant times the speed of light. For the conversion between energy and photons.
        
        """
        self.const = StructuredConfiguration()
        self.const.setFile("const","SED.const.config.yaml")
        self.const["hc"] = 1.98644521e-8 # erg angstrom
        self.Caches["CONST"] = ConfigCache(self.const,filename=self.config["Caches.const"])
        
    
    
    @description("Setting up lenslets")
    @depends("setup-config")
    def setup_lenslets(self):
       """This function loads data about lenslet positions, and thier dispersion through the prism. The data are original produced by ZEMAX. This function reads the Zeemax data directly and then cleans the data in certain ways, preparing it for use later in the system.
       
       **Command Name:** ``*setup-lenslets``
       
       Cleaning actions taken:
       
       - Indexes (ix) become integers
       - Wavelenghts are converted to SI units (m) instead of microns
       - Center of the lenslet array is calculated
       - Lenslets are validated. See :meth:`SEDMachine.Lenslet.Lenslet.valid`.
        
       **Data expected from ZEMAX**:
       
       The data are given as a list of spots. Each lenslet will have many spots. The variables below are listed in order.
       
       :var ix: index of the lenslet
       :var xps: Pupil position in the x-direction in mm
       :var yps: Pupil position in the y-direction in mm
       :var lams: Wavelength of this spot in microns
       :var xcs: Camera image position of the spot in the x direction in mm.
       :var ycs: Camera image position of the spot in the y direction in mm.
       :var xls: Camera image position of the next R=100 resolution element in the x direction in mm.
       :var yls: Camera image position of the next R=100 resolution element in the y direction in mm.
       :var xas: Camera image position of the major axis extent of the telescope image in the x direction in mm.
       :var yas: Camera image position of the major axis extent of the telescope image in the y direction in mm.
       :var xbs: Camera image position of the minor axis extent of the telescope image in the x direction in mm.
       :var ybs: Camera image position of the minor axis extent of the telescope image in the y direction in mm.
       :var rs: Instantaneous resolution of this position.
       
      
       """
       # Load Lenslet Specification File
       self.log.debug("Opening filename %s" % self.config["Instrument.files.lenslets"])
       ix, xps, yps, lams, xcs, ycs, xls, yls, xas, yas, xbs, ybs, rs = np.genfromtxt(self.dir_filename("Data",self.config["Instrument.files.lenslets"]),skip_header=1,comments="#",unpack=True)
       # This data describes the following:
       # ix - Index (number)
       # xps - Pupil position in the x-direction
       # yps - Pupil position in the y-direction
       # lams - wavelengths for this position
       # xs - X position (in mm, offest from top right corner) of this wavelength
       # ys - Y Position (in mm, offset from top right corner) of this wavelength
        
       # Correctly Type Lenslet Specification Data
       ix = ix.astype(np.int) #Indicies should always be integers
        
       lams *= 1e-6 # Convert wavelength to SI units (m)
        
       # This simply generates a list of all of the lenslets
       self.lensletIndex = np.unique(ix)
        
       # Determine the center of the whole system by finding the x position that is closest to 0,0 in pupil position
       cntix = np.argmin(xps**2 + yps**2)
       self.center = ((xcs[cntix] + (self.config["Instrument.image.size.mm"]/2))* self.config["Instrument.convert.mmtopx"], (ycs[cntix] + (self.config["Instrument.image.size.mm"]/2)) * self.config["Instrument.convert.mmtopx"])
       
        
       # Progress bar for lenslet creation and validation
      
       total = len(self.lensletIndex)
       self._start_progress_bar(total,"green")
        
       # Variables for lenslet use
       FileName = "%(Partials)s/%(name)s%(ext)s" % dict(name="Lenslets-raw",ext=".dat",**self.config["Dirs"])
       with open(FileName,'w') as stream:
           for idx in self.lensletIndex:
               select = idx == ix
               lenslet = Lenslet(xps[select],yps[select],lams[select],idx,xcs[select], ycs[select],xls[select], yls[select],  xas[select], yas[select], xbs[select], ybs[select], rs[select],self.config,self.Caches)
               if lenslet.valid(strict=self.config["Instrument.Lenslets.strict"]):
                   self.lenslets[idx] = lenslet
                   stream.write(lenslet.introspect())
               self.progress += 1
               self.progressbar.update(self.progress)
       self.lensletIndex = np.asarray(self.lenslets.keys())
       self._end_progress_bar()
       
       
       # Central Lenslet Output
       FileName = "%(Partials)s/%(name)s%(ext)s" % dict(name="center-raw",ext=".dat",**self.config["Dirs"])
       with open(FileName,'w') as stream:
           lenslet = self.lenslets[ix[cntix]]
           stream.write(lenslet.introspect())
       
       
       if "start" in self.config["Lenslets"]:
           self.lensletIndex = self.lensletIndex[self.config["Lenslets.start"]:]
       if "number" in self.config["Lenslets"]:
           self.lensletIndex = self.lensletIndex[:self.config["Lenslets.number"]]
       if "radius" in self.config["Lenslets"] and "position" in self.config["Lenslets"]:
           xp,yp = self.config["Lenslets.position.x"],self.config["Lenslets.position.y"]
           xps,yps = np.array([l.ps[0] for l in self.lenslets.values()]).T
           distances = np.sqrt((xps-xp)**2.0 + (yps-yp)**2.0)
           include = distances <= self.config["Lenslets.radius"]
           self.lensletIndex = self.lensletIndex[include]
       self.total = len(self.lensletIndex)
       self.lenslets = {x:self.lenslets[x] for x in self.lensletIndex}
       for ix,lx in enumerate(self.lenslets.values()):
           lx.idx = ix
    
    
    @description("Creating blank frame")
    @depends("setup-config")
    def setup_blank(self):
        """Establish a blank image of zeros in every position. The image is established as ``np.int32`` values.
        
        **Command Name:** ``*setup-blank``
        """
        self["Blank"] = np.zeros((self.config["Instrument.image.size.px"],self.config["Instrument.image.size.px"])).astype(np.int32)
        
    
    @description("Creating blank frame with a single one")
    @depends("setup-config")
    @replaces("setup-blank")
    def setup_dummy_blank(self):
        """Setup a dummy blank image with a single value of 1.0 in the center of the image. This is useful for a sanity-check of the scattered light calculation.
        
        **Command Name:** ``*setup-blank-d``
        """
        blank = np.zeros((self.config["Instrument.image.size.px"],self.config["Instrument.image.size.px"])).astype(np.int32)
        center = np.int(self.config["Instrument.image.size.px"]/2.0)
        blank[center,center] = 1.0
        self["Blank"] = blank
        
    def setup_source(self):
        """Set up a source for use with this system.
        
        **Command Name:** ``*setup-source``
        
        .. Warning::
            This method is not ready for use yet. It requires some concept of the wavelength data for a data-cube. Extracting the wavelength calibration may be non-trivial."""
        
        self.log.warning("Stage 'setup-source' not ready yet, doing nothing!")
        return
        
        Source = ImageStack()
        Source.read(self.dir_filename("Data",self.config["Source.CubeName"]))
        data = Source.data()
        shape = data.shape
        
        # Progress bar for lenslet creation and validation
        PBar = arpytools.progressbar.ProgressBar(color="green")
        finished = 0.0
        total = shape[0] * shape[1]
        self.log.useConsole(False)
        PBar.render(0,"L:%4s,%4s %6d/%-6d" % ("","",finished,total))
        
        self.SourcePixels = []
        for i in range(shape[0]):
            for j in range(shape[1]):
                progress = int((finished/float(total)) * 100)
                finished += 1
                self.SourcePixels.append(SourcePixel(i,j,data=data[i,j],label="Source Pixel %d,%d" % (i,j),config=self.config,num=finished))
                PBar.render(progress,"L:%4d,%4d %6d/%-6d" % (i,j,finished,total))
        
        PBar.render(1,"L:%4s,%4s %6d/%-6d" % ("Done","",finished,total))
        
    
    
    @description("Creating a simple source object")
    @depends("setup-config","setup-constants")
    @replaces("setup-source")
    def setup_simple_source(self):
        """Creates a simple-source object. The simple source is placed in three fixed positions on the image plane. Each source is a uniform distribution of the provided spectrum.
        
        **Command Name:** ``*simple-source``        
        
        The spectrum is set up using the source configuration values::
            
            Source:
              Filename: Data/SNIa.R1000.dat
        
        There is no amplification applied to the source. Sources are expected to be in cgs units during input.
            
        """        
        WL,FL = np.genfromtxt(self.dir_filename("Data",self.config["Source.Filename"]),unpack=True,comments="#")
        FL /= self.const["hc"] / WL
        FL *= 1e10 #Spectrum was per Angstrom, should now be per Meter
        WL *= 1e-10
        self.spectra.save(InterpolatedSpectrum(np.array([WL,FL]),self.config["Source.Filename"],method="resolve_and_integrate"))
        self.spectra.save(InterpolatedSpectrum(np.array([WL,FL]),self.config["Source.Filename"]+" (O)",method="resolve_and_integrate"),select=False)
        self.spectra.save(SpectraFrame(np.array([WL,FL]),"Original"),select=False)
        self.SourcePixels = [SourcePixel(-0.13,0,data=np.array([WL,FL]),label="Source Pixel",config=self.config,num=1),SourcePixel(0,0,data=np.array([WL,FL]),label="Source Pixel",config=self.config,num=2),SourcePixel(0.05,0.05,data=np.array([WL,FL]),label="Source Pixel",config=self.config,num=3)]
        for ix,px in enumerate(self.SourcePixels):
            px.idx = ix
    
    
    @description("Setting up Cameras")
    @depends("setup-config")
    def setup_cameras(self):
        """Set up camera configuration values. Camera configuratiosn can be slightly dynamic, so some values are copied from others.
        
        .. Note:: 
            Camera configurations do not contain QE values (unlike the SEDSpec simulators).
        
        """
        self.config["Instrument.Cameras.PI-fast"] = self.config["Instrument.Cameras.PI"]
        self.config["Instrument.Cameras.PI-fast.RN"] = 12
        self.config["Instrument.Cameras.PI-fast.readtime"] = 2.265

        self.config["Instrument.Cameras.Andor-fast"] = self.config["Instrument.Cameras.Andor"]
        self.config["Instrument.Cameras.Andor-fast.RN"] = 11.7
        self.config["Instrument.Cameras.Andor-fast.readtime"] = 1.398
        
    
    
    @description("Setting up Sky Sources")
    @depends("setup-config","setup-constants")
    def setup_sky(self):
        """Setup sky spectrum from a file. Sky spectrum setup also includes moon spectrum, atmospheric extinction setup and quantum efficiency setup for cameras.
        
        **Command Name:** ``*setup-sky``
        
        Configuration for this stage is included in both the ``Instrument`` values and ``Observation``::
            
            Observation:
              Background:
                Atmosphere: Atmosph
                Files:
                  Atmosph:
                    Amplifier: 1
                    Filename: Data/atmosphere.fits
                  PalSky:
                    Amplifier: 1
                    Filename: Data/PalSky.fits
                Sky: PalSky
              Moon:
                Phase: 0.45
              airmass: 1
            Instrument:
              Thpt:
                File: Data/thpt.npy
                Type: prism_pi
            
            
        These configuration values allow for multiple possibilities for sky and atmosphere spectra. See the default configuration values found with::
            
            $ SEDMsim --dump *none
            
        for more information.
            
        """
        
        
        # Sky Data (From sim_pdr.py by Nick, regenerated using SEDTools module's make_files.py script)
        # Each sky spectrum is saved in a FITS file for easy recall as a spectrum object.
        self.SKYData = SpectraStack()
        for label,d in self.config["Observation.Background.Files"].iteritems():
            self.SKYData.load(self.dir_filename("Data",d["Filename"]),framename=label)
                
        
        # Moon phase adjustments. These moon phase attenuation values are for different filter bands.
        # The intermediate wavelengths are accounted for using a polyfit
        # the result is a function which takes phase and wavelength, and outputs an attenuation...
        
        # See derivation on pg 83 of SED NB 1 (20 July 2011)
        moon_phase = np.array([0., 0.08, 0.16, 0.24, 0.32, 0.40, 0.50])
        moon_g = np.array([2e-17, 2.1e-17, 2.15e-17, 2.3e-17, 5.3e-17, 1.7e-16, 3.2e-16])
        moon_r = np.array([2.3e-17,2.3e-17,2.3e-17,3.3e-17,3.5e-17,8.3e-17,1.3e-16])
        moon_i = np.array([2.8e-17,3.0e-17,3.0e-17,3.3e-17,3.8e-17,7.0e-17,9.0e-17])
        
        sky_ls = np.array([4868., 6290., 7706., 10000]) * 1e-10
        
        gm = moon_g - moon_g[0]
        rm = moon_r - moon_r[0]
        im = moon_i - moon_i[0]
        zm = im
        
        fluxes = np.array([gm, rm, im, zm])
        
        moon_specs = [ scipy.interpolate.interp1d(moon_phase,fl) for fl in fluxes]            
        mfl = []
        mls = []
        for i in xrange(len(moon_specs)):
            mfl.append(moon_specs[i](self.config["Observation.Moon.Phase"]))
            mls.append(sky_ls[i])
            
        
        
        # Throughputs are generated from Nick's simulation scripts in throughput.py
        # They are simply re-read here.
        thpts = np.load(self.dir_filename("Data",self.config["Instrument.Thpt.File"]))[0]
        WL = thpts["lambda"]* 1e-10
        self.qe["prism_pi"] = InterpolatedSpectrum(np.array([WL, thpts["thpt-prism-PI"]]),"PI Prism")
        self.qe["prism_andor"] = InterpolatedSpectrum(np.array([WL, thpts["thpt-prism-Andor"]]),"Andor Prism")
        self.qe["grating"] = InterpolatedSpectrum(np.array([WL, thpts["thpt-grating"]]),"Grating")
        self.qe.select(self.config["Instrument.Thpt.Type"])
        
        # Set up extinction and airmass term.
        WL,EX = self.SKYData.data(self.config["Observation.Background.Atmosphere"])
        FL = 10**(-EX*self.config["Observation.airmass"]/2.5)
        WL *= 1e-10
        self.sky.save(InterpolatedSpectrum(np.array([WL,FL]),"Atmosphere"),select=False)
        
        
        # This calculation fixes the units of the TurnroseSKY values
        # I'm not sure what these units are doing, but we will leave them here for now.
        WL,FL = self.SKYData.data(self.config["Observation.Background.Sky"])
        FL *=  self.config["Observation.Background.Files"][self.config["Observation.Background.Sky"]]["Amplifier"]
        FL /= self.const["hc"] / WL
        FL *= 1e10 #Spectrum was per Angstrom, should now be per Meter
        
        
        M_FL = InterpolatedSpectrum(np.array([mls,mfl]),"Moon Phase %s" % i,method='polyfit')(wavelengths=WL*1e-10)[1]
        M_FL /= self.const["hc"] / WL
        M_FL *= 1e10
        
        WL *= 1e-10
        
        
        self.sky.save(InterpolatedSpectrum(np.array([WL,M_FL]),"Moon",method='resolve_and_integrate'),select=False)
        self.sky.save(InterpolatedSpectrum(np.array([WL,FL]),"SkySpectrum (O)",method="resolve_and_integrate"),select=False)
        FL += M_FL
        self.sky.save(InterpolatedSpectrum(np.array([WL,FL]),"SkySpectrum",method="resolve_and_integrate"))
        
    
    @description("Applying atmospheric extinction term")
    @depends("setup-sky")
    def apply_atmosphere(self):
        """Apply the atmospheric extinction term."""
        self.map_over_lenslets(self._apply_atmosphere,color=False)
        self.spectra.save( self.spectra.frame() * self.sky.frame("Atmosphere") )
        
        
    def _apply_atmosphere(self,lenslet):
        """docstring for _apply_atmosphere"""
        lenslet.spectrum *= self.sky.frame("Atmosphere")
        
    
    @description("Setting up lenslet hexagons")
    @depends("setup-lenslets")
    def setup_hexagons(self):
        """Make the lenslet hexagons"""
        self.map_over_lenslets(lambda l: l.make_hexagon(),color="green")
        
    
    @description("Making source pixels")
    @depends("setup-source")
    def setup_source_pixels(self):
        """Setup source pixels"""
        self.map_over_pixels(lambda p: p.make_pixel_square(),color="green")
       
    
    @description("Applying Sky Spectrum to Source")
    @depends("setup-sky")
    def apply_sky(self):
        """Apply Sky Spectrum to each lenslet"""
        self.map_over_lenslets(self._apply_sky_spectrum,color=False)
        self.spectra.save(self.spectra.frame() + self.sky.frame())
        
        
    def _apply_sky_spectrum(self,lenslet):
        """docstring for _apply_sky_spectrum"""
        lenslet.spectrum += self.sky.frame()
        
    def _apply_qe_spectrum(self,lenslet):
        """Apply qe to each lenslet"""
        lenslet.spectrum *= self.qe.frame() * self.config["Instrument.Tel.area"] * self.config["Observation.exposure"] * self.config["Instrument.eADU"]
        
    
    @description("Applying Quantum Efficiency Functions")
    @depends("setup-sky")
    def apply_qe(self):
        """Apply the instrument quantum efficiency"""
        self.map_over_lenslets(self._apply_qe_spectrum,color=False)
        self.sky.save(self.sky.frame() * self.config["Instrument.Tel.area"] * self.config["Observation.exposure"] * self.config["Instrument.eADU"],"Sky Mul")
        self.sky.save(self.sky.frame() * self.qe.frame())
        self.spectra.save(self.spectra.frame() * self.config["Instrument.Tel.area"] * self.config["Observation.exposure"] * self.config["Instrument.eADU"],"Spec Mul")
        self.spectra.save(self.spectra.frame() * self.qe.frame())
        self.sky.save(self.sky.frame("Moon") * self.config["Instrument.Tel.area"] * self.config["Observation.exposure"] * self.config["Instrument.eADU"],"Moon Mul",select=False)
        self.sky.save(self.sky.frame("Moon Mul") * self.qe.frame(),"Moon QE",select=False)
        
    
    @description("Performing geometric resample")
    @depends("setup-source-pixels","setup-hexagons")
    def geometric_resample(self):
        """docstring for fname"""
        n = len(self.SourcePixels)
        m = len(self.lenslets)
        self.map_over_lenslets(lambda l:l.setup_crosstalk(n),color=False)
        self.map_over_pixels(lambda p:p.setup_crosstalk(m),color=False)
        self.map_over_lenslets(lambda l:self.map_over_pixels(lambda p:l.find_crosstalk(p),color=False),color="green")
        
    
    @description("Setting up calibration source")
    @depends("setup-config","setup-constants")
    def setup_line_list(self):
        """Set up a line-list based spectrum for wavelength calibration."""
        linelist = np.asarray(np.genfromtxt(self.dir_filename("Data",self.config["Source.Lines.List"]),comments="#"))
        if linelist.ndim < 1:
            linelist = linelist.flatten()
        CalSpec = FlatSpectrum(0.0)
        sigma = self.config["Source.Lines.sigma"]
        for line in linelist:
            if self.config["Source.Lines.peaks"]:
                value = line[1]
                center = line[0]
            else:
                value = self.config["Source.Lines.value"]
                center = line
            CalSpec += GaussianSpectrum(center,sigma,value,"Line %g" % line)
        self.spectra.save(UnitarySpectrum(CalSpec,method='resolve_and_integrate',label="Calibration Lamp"),select=False)
        
    
    @description("Using calibration lamp source")
    @help("Use a calibration lamp source")
    @depends("setup-lenslets","setup-lines","setup")
    @replaces("setup-source","geometric-resample","setup-source-pixels","apply-sky","apply-atmosphere")
    def line_source(self):
        """Use the line spectrum only"""
        self.config["Output.Label"] += "-cal-"
        self.replace_source(self.spectra.frame("Calibration Lamp"))
        
    
    @description("Using only Sky spectrum")
    @description("Use only Sky spectrum")
    @depends("setup-sky","apply-sky","apply-qe","apply-atmosphere","setup-lenslets","setup")
    @replaces("setup-source","geometric-resample","setup-source-pixels")
    def sky_source(self):
        """Use the sky spectrum only"""
        self.config["Output.Label"] += "-sky-"
        self.replace_source(self.sky.frame())

    
    @ignore
    def replace_source(self,spectrum):
        """Replace the default file-source with a flat spectrum"""
        self.spectra.save(spectrum,clobber=True)
        self.map_over_lenslets(self._replace_source,color=False)
    
    
    @ignore
    def _replace_source(self,lenslet):
        """docstring for _flat_source"""
        lenslet.spectrum = self.spectra.frame()
            
    
    
    
    @description("Using flat source")
    @help("Make a constant value source")
    @depends("setup-lenslets","setup")
    @replaces("setup-source","geometric-resample","setup-source-pixels","apply-sky","apply-atmosphere")
    def flat_source(self):
        """Replace the default file-source with a flat spectrum"""
        self.config["Output.Label"] += "-flat-"
        self.replace_source(FlatSpectrum(self.config["Source.Flat.value"]))
        
    
    @description("Calculating dispersion for each lenslet")
    @help("Calculate lenslet dispersion")
    @depends("setup-lenslets","setup-caches")
    def lenslet_dispersion(self):
        """Calculate the dispersion for each lenslet"""
        self.map_over_lenslets(lambda l: l.find_dispersion(),color="blue")
        
    
    @description("Tracing lenslet spectra dispersion")
    @help("Trace lenslet spectra dispersion")
    @depends("dispersion","setup-caches","setup-lenslets","apply-sky","apply-qe","apply-atmosphere")
    def lenslet_trace(self):
        """Trace out each lenslet"""
        self.map_over_lenslets(lambda l: l.get_trace(l.spectrum),color="blue")

    
    @include
    @description("Placing lenslet spectra")
    @help("Place subimages")
    @depends("trace")
    def lenslet_place(self):
        """Place each spectrum into the subimage"""
        self.map_over_lenslets(self._lenslet_place,color="yellow")
        
    
    @ignore
    def _lenslet_place(self,l):
        """docstring for _lenslet_place"""
        l.place_trace(self.get_conv)
        l.write_subimage()
    
    
    @include
    @description("Merging subimages")
    @depends("setup-blank","setup-lenslets")
    def image_merge(self):
        """Merge subimages into master image"""
        self.select("Blank")
        self.save(self.frame(),"Merge")
        self.map_over_lenslets(self._lenslet_merge,color="yellow")
        self.select("Merge")
        
    
    
    @ignore
    def _lenslet_merge(self,lenslet):
        """Merge a single lenslet into the master image"""
        lenslet.read_subimage()
        lenslet.bin_subimage()
        self.place(lenslet.data(),lenslet.subcorner)
        lenslet.clear(delete=True)
    
    
    
    @description("Setting up scattered light calculations")
    @depends("setup-config","setup-blank")
    def setup_scatter(self):
        """Sets up scattered light level"""
        
        area = np.zeros((self.config["Instrument.ccd.size.px"]+1,self.config["Instrument.ccd.size.px"]+1))
        
        size = self.config["Instrument.ccd.size.px"]/2
        for v in self.config["Instrument.Scatter.Kernels"].values():
            if v["type"] == "Gaussian":
                n = v["mag"] * self.gauss_kern(v["stdev"],size,enlarge=False,normalize=False)
                self.log.debug(npArrayInfo(n,"Gauss Kernel"))
            area += n
        
        state = self.framename
        self.save(area,"Scatter")        
        self.select(state)
        
    
    
    @include
    @help("Crop Final Image")
    @description("Cropping image to CCD size")
    @depends("setup-blank","setup-lenslets")
    def ccd_crop(self):
        """Crops the image to the appropriate ccd size"""
        x,y = self.center
        size = self.config["Instrument.ccd.size.px"] / 2.0
        self.crop(x,y,size)
        
    
    @description("Setting up Dark/Bias frames")
    @depends("setup-config","setup-cameras")
    def setup_noise(self):
        """Makes noise masks"""
        
        state = self.framename
        
        read_noise = self.config["Instrument.Cameras"][self.config["Instrument.Cameras.Selected"]]["RN"]
        dark_noise = self.config["Instrument.Cameras"][self.config["Instrument.Cameras.Selected"]]["DC"] * self.config["Observation.exposure"]
        
        self.generate_poisson_noise("Read",read_noise,self.config["Observation.number"])
        
        self.generate_poisson_noise("Dark",dark_noise)
        
        self.select(state)
    
        
    
    @include
    @help("Add Dark/Bias noise to image")
    @description("Adding Dark/Bias noise")
    @depends("crop","setup-noise")
    def apply_noise(self):
        """Apply the noise masks to the target image label"""
        
        dark = self.data("Dark")
        bias = self.data("Read")
        
        data = self.data()
        
        data += dark + bias
        
        self.save(data,"Noisy",clobber=True)
    
    
    @include
    @help("Add scattered light noise to image")
    @description("Adding scatter noise")
    @depends("crop","setup-scatter")
    def apply_scatter(self):
        """Apply the scattered light frame"""
        
        
        scatter = self.data("Scatter")
        
        self.log.debug(npArrayInfo(scatter,"Scatter"))
        
        
        data = self.data()
        self.log.debug(npArrayInfo(data,"Data \'%s\'" % self.framename))
        
        
        if self.config["Instrument.Scatter.FFT"]:
            result = sp.signal.fftconvolve(data,scatter,mode='same')
        else:
            result = sp.signal.convolve(data,scatter,mode='same')
        
        result *= self.config["Instrument.Scatter.Amplifier"]
        
        self.log.debug(npArrayInfo(result,"Scattered Light"))
        end = data + result[:-1,:-1]
        
        self.save(result,"ScatterOnly")
        self.save(end,"Scattered",clobber=True)
        
    
    
    @include
    @help("Transpose the image")
    @description("Transposing Image")
    @depends("crop")
    def transpose(self):
        """transpose the final image"""
        data = self.data()

        self.save(data.T.astype(np.int16),"Transposed",clobber=True)
        
    
    
    
    
    
    @include
    @help("Save image to file")
    @description("Saving image to disk")
    def save_file(self):
        """Saves the file"""
        self.Filename = "%(Output)s/%(label)s-%(date)s.%(fmt)s" % dict(label=self.config["Output.Label"],date=time.strftime("%Y-%m-%d"), fmt=self.config["Output.Format"], **self.config["Dirs"] )
        self.write(self.Filename,frames=[self.framename],clobber=True)
        self.log.info("Wrote %s" % self.Filename)
        self.Filename = "%(Output)s/%(label)s-deep-%(date)s.%(fmt)s" % dict(label=self.config["Output.Label"],date=time.strftime("%Y-%m-%d"), fmt=self.config["Output.Format"], **self.config["Dirs"] )
        self.write(self.Filename,clobber=True)
        self.log.info("Wrote %s" % self.Filename)
    
    
    ################################
    ## IMAGE MANAGEMENT FUNCTIONS ##
    ################################
    
    @ignore
    def place(self,img,corner):
        """Place the given AstroObject.AnalyticSpectra.AnalyticSpectrum onto the SEDMachine Image"""
        
        xstart = corner[0]
        xend = xstart + img.shape[0]
        ystart = corner[1]
        yend = ystart + img.shape[1]
        data = self.data()
        
        if data.shape[0] < xend or data.shape[1] < yend:
            raise SEDLimits
        
        if xstart < 0 or ystart < 0:
            raise SEDLimits
        
        if xend < 0 or yend < 0:
            raise SEDLimits
        
        data[xstart:xend,ystart:yend] += img
        self.save(data,self.framename,clobber=True)    
        
    
    
    @ignore
    def crop(self,x,y,xsize,ysize=None,label=None):
        """Crops the provided image to twice the specified size, centered around the x and y coordinates provided."""
        if not ysize:
            ysize = xsize
        cropped = self.d[x-xsize:x+xsize,y-ysize:y+ysize]
        self.log.debug("Cropped and Saved Image")
        if label == None:
            label = "Cropped"
        self.save(cropped,label,clobber=True)
    
    
    #######################
    ## DEBUGGING METHODS ##
    #######################
    
    @depends("dispersion")
    @description("Plotting dispersion for each lenslet")
    def plot_dispersion_data(self):
        """Outputs dispersion debugging data"""
        self.map_over_lenslets(lambda l: l.plot_dispersion(),color="cyan")
    
    @description("Plotting trace data for each lenslet")
    @depends("trace")
    def plot_trace_data(self):
        """Outputs plots about each lenslet trace"""
        self.map_over_lenslets(lambda l: l.plot_trace(),color="cyan")
        
    def plot_spectrum_data(self):
        """Outputs plots about each lenslet trace"""
        self.map_over_lenslets(lambda l: l.plot_spectrum(),color="cyan")
    
    @description("Plotting lenslet positions")
    @depends("setup-lenslets")
    def plot_lenslet_data(self):
        """Outputs the lenslet data"""
        plt.figure()
        plt.clf()
        self.log.info("Plotting lenslet arc positions in Camera (x,y) space")
        FileName = "%(Partials)s/Lenslet-cxy%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        self.map_over_lenslets(lambda l: plt.plot(l.xcs,l.ycs,linestyle='-'),color="cyan")
        plt.title("Lenslet c-xy positions")
        plt.savefig(FileName)
        
        plt.clf()
        self.log.info("Plotting lenslet Pupil positions in Pupil (x,y) space")
        FileName = "%(Partials)s/Lenslet-pxy%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        self.map_over_lenslets(lambda l: plt.plot(l.ps.T[0],l.ps.T[1],marker='.'),color="cyan")
            
            
        plt.title("Lenslet p-xy positions")
        plt.savefig(FileName)
        
        plt.clf()
        
    
    def plot_lenslet_raw(self):
        """Plot lesnlet raw data"""
        self.map_over_lenslets(lambda l: l.plot_raw_data(),color="cyan")
    
    @description("Plotting Spectrum Tests")
    @depends("setup-sky","setup-source","apply-sky","apply-qe","apply-atmosphere")
    def compare_methods(self):
        """A plot for comparing resolving methods"""
        WL = self.config["Instrument.wavelengths.values"]
        RS = self.config["Instrument.wavelengths.resolutions"]
        DWL,DRS = self.get_resolution_spectrum(np.min(WL),np.max(WL),1000)
        SWL,SRS = self.get_resolution_spectrum(np.min(WL),np.max(WL),500)
        plt.figure()
        plt.title("Resolutions")
        plt.plot(WL*1e6,RS,'g.')
        plt.plot(DWL*1e6,DRS,'r.')
        plt.plot(SWL*1e6,SRS,'b.')
        plt.axis(expandLim(plt.axis()))
        FileName = "%(Partials)s/Test-Spectrum-Res%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)

        plt.clf()
        plt.title("Resolve \& Resample Tests")
        self.log.debug(npArrayInfo(WL,"Wavelength for Sky Plot"))
        
        Results = {}
        master = "H Iq"
        
        for wl,rs,Label in zip([WL,DWL,SWL],[RS,DRS,SRS],["L","M","H"]):
            identity = "%s Iq" % Label                
            wl,FL = self.spectra.frame(self.config["Source.Filename"]+" (O)")(wavelengths=wl,resolution=rs,method="integrate_quad")
            Results[identity] = np.sum(FL)
            plt.semilogy(wl*1e6,FL,'g.',linestyle="-",label="%s $\Sigma$%.2e" % (identity,Results[identity]))
            
            identity = "%s Ih" % Label                
            wl,FL = self.spectra.frame(self.config["Source.Filename"]+" (O)")(wavelengths=wl,resolution=rs,method="integrate_hist")
            Results[identity] = np.sum(FL)
            plt.semilogy(wl*1e6,FL,'b.',linestyle="-",label="%s $\Sigma$%.2e" % (identity,Results[identity]))
            
            identity = "%s RR" % Label                
            wl,FL = self.spectra.frame(self.config["Source.Filename"]+" (O)")(wavelengths=wl,resolution=rs,method="resolve_and_integrate")
            Results[identity] = np.sum(FL)
            plt.semilogy(wl*1e6,FL,'r.',linestyle="-",label="%s $\Sigma$%.2e" % (identity,Results[identity]))
        
        text = "Error relative to %s in integration and resolution spectrum methods:\n" % master
        mv = Results[master]
        for Label in ["L","M","H"]:
            for method in ["Iq","Ih","RR"]:
                identity = "%(label)s %(method)s" % {'label':Label,'method':method}
                val = Results[identity]
                perr = (np.abs(val - mv) / mv) * 100.0
                line = "%(identity)s $\Sigma = $%(value).3e Error %(perr).2f \%%\n" % {'identity':identity,'value':val,'perr':perr}
                self.log.debug(line)
                if identity == master:
                    line = line.rstrip("\n")
                    line += " MASTER\n"
                text += line
        text = text.rstrip("\n")
        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("Flux (Photons)")
        plt.legend(loc=3, mode="expand", borderaxespad=0.,ncol=3)
        FileName = "%(Partials)s/Test-R-and-R%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
        
        
        
        plt.title("Resolve \& Resample Test Results")
        self.log.debug(text)
        ax = plt.gca()
        plt.text(0.2, 0.5,text,
            horizontalalignment='left',
            verticalalignment='center',
            transform = ax.transAxes)
        FileName = "%(Partials)s/Test-R-and-R-vals%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        
        
        plt.clf()
        plt.title("Resample Tests")
        self.log.debug(npArrayInfo(WL,"Wavelength for Sky Plot"))
        
        WL,FL = self.spectra.frame(self.config["Source.Filename"]+" (O)")(wavelengths=WL,method="interpolate")
        plt.semilogy(WL*1e6,FL,'-',label="Interpolate")
        
        SWL,FL = self.spectra.frame(self.config["Source.Filename"]+" (O)")(wavelengths=SWL,resolution=SRS,method="resample")
        plt.semilogy(SWL*1e6,FL,'-',label="SR Resample")
        
        DWL,FL = self.spectra.frame(self.config["Source.Filename"]+" (O)")(wavelengths=DWL,resolution=DRS,method="resample")
        plt.semilogy(DWL*1e6,FL,'-',label="HR Resample")
        
        WL,FL = self.spectra.frame(self.config["Source.Filename"]+" (O)")(wavelengths=WL,resolution=RS,method="resample")
        plt.semilogy(WL*1e6,FL,'-',label="LR Resample")
        
        


        
        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("Flux (Photons)")
        plt.legend(loc=4)
        FileName = "%(Partials)s/Test-Resample%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
        
        
    @description("Plotting Sky Spectrum")
    @depends("setup-sky","apply-sky","apply-qe","plot-sky-o")
    def plot_sky(self):
        """Plot sky spectrum"""
        WL = self.config["Instrument.wavelengths.values"]
        RS = self.config["Instrument.wavelengths.resolutions"]
        plt.figure()
        plt.title("Resolution")
        plt.plot(WL*1e6,RS,'g.')
        FileName = "%(Partials)s/Sky-Spectrum-Res%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)

        plt.clf()
        plt.title("Sky Spectrum")

        WL,FL = self.sky.frame()(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'b.',linestyle='-',label="Sky + Moon (qe)")
        
        WL,FL = self.sky.frame("Sky Mul")(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'m.',linestyle='-',label="Sky + Moon")
        
        WL,FL = (self.sky.frame("Sky Mul") - self.sky.frame("Moon Mul"))(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'g.',linestyle='-',label="Sky")
        
        axis = plt.axis()
        
        WL,FL = self.sky.frame("Moon Mul")(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'y.',linestyle='-',label="Moon",zorder=0.5)
        
        WL,FL = self.sky.frame("Moon QE")(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'c.',linestyle='-',label="Moon (qe)",zorder=0.5)
        
        # plt.axis(axis)
        plt.legend(loc=2)

        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("Flux (Photons)")
        FileName = "%(Partials)s/Sky-Spectrum%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
    
    @description("Plotting Original Sky Spectrum")
    @depends("setup-sky")
    def plot_sky_original(self):
        """Plot sky spectrum"""
        WL = self.config["Instrument.wavelengths.values"]
        RS = self.config["Instrument.wavelengths.resolutions"]
        plt.figure()
        plt.title("Resolution")
        plt.plot(WL*1e6,RS,'g.')
        FileName = "%(Partials)s/Sky-Spectrum-Res%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)

        plt.clf()
        plt.title("Sky Spectrum")

        WL,FL = self.sky.frame("SkySpectrum (O)")(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'g.',linestyle='-',label="Sky")
        
        plt.legend(loc=2)

        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("Flux (Photons)")
        FileName = "%(Partials)s/Sky-Original-Spectrum%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
    
    
    @description("Plotting Original Source Spectrum")
    @depends("setup-lines")
    def plot_original_calibration(self):
        """Plot the original calibration source"""
        WL = self.config["Instrument.wavelengths.values"]
        RS = self.config["Instrument.wavelengths.resolutions"]
        plt.figure()
        plt.title("Calibration Spectrum")
        WL,FL = self.spectra.frame("Calibration Source")(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'r.',linestyle='-',label="Source")
        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("Flux (Photons)")
        plt.legend(loc=4)
        FileName = "%(Partials)s/Cal-Original-Spectrum%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
        
    
    def plot_original_source(self):
        """Plot the original source spectrum only"""
        WL = self.config["Instrument.wavelengths.values"]
        RS = self.config["Instrument.wavelengths.resolutions"]
        plt.figure()
        plt.title("Resolution")
        plt.plot(WL*1e6,RS,'g.')
        FileName = "%(Partials)s/Source-Spectrum-Res%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)

        plt.clf()
        plt.title("Source Spectrum")
        WL,FL = self.spectra.frame(self.config["Source.Filename"]+" (O)")(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'r.',linestyle='-',label="Source")
        
        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("Flux (Photons)")
        plt.legend(loc=4)
        FileName = "%(Partials)s/Source-Original-Spectrum%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
        
    
    @description("Plotting Source Spectrum")
    @depends("setup-sky","setup-source","apply-sky","apply-qe","apply-atmosphere","plot-source-o")
    def plot_source(self):
        """Plot the source spectrum"""
        WL = self.config["Instrument.wavelengths.values"]
        RS = self.config["Instrument.wavelengths.resolutions"]
        plt.figure()
        plt.title("Resolution")
        plt.plot(WL*1e6,RS,'g.',label="Requested R")
        GWL,GFL = self.spectra.data("Original")
        plt.plot(GWL[:-1]*1e6,GWL[:-1]/np.diff(GWL),label="Given R")
        plt.legend()
        FileName = "%(Partials)s/Source-Spectrum-Res%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)

        plt.clf()
        plt.title("Source Spectrum")
        self.log.debug(npArrayInfo(WL,"Wavelength for Sky Plot"))
        
        WL,FL = self.spectra.frame()(wavelengths=WL,resolution=RS)
        self.log.debug(npArrayInfo(WL,"Wavelength from Source Plot"))
        self.log.debug(npArrayInfo(FL,"Flux from Source Plot"))
        plt.semilogy(WL*1e6,FL,'b.',linestyle='-',label="Combined")
        
        WL,FL = (self.spectra.frame() / self.sky.frame("Atmosphere"))(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'r.',linestyle='-',label="Source (qe)")
        
        WL,FL = self.sky.frame()(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'g-',linestyle='-',label="Sky (qe)")
        
        WL,FL = self.sky.frame("Sky Mul")(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'g-',linestyle='-',label="Sky + Moon",zorder=0.5)
        
        
        WL,FL = (self.sky.frame("Sky Mul") - self.sky.frame("Moon Mul"))(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'y-',linestyle='-',label="Sky")
        
        WL,FL = (self.spectra.frame("Spec Mul") - self.sky.frame("Moon Mul"))(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'m.',linestyle='-',label="Source")
        
        axis = plt.axis()
        
        
        WL,FL = self.sky.frame("Moon Mul")(wavelengths=WL,resolution=RS)
        self.log.debug(npArrayInfo(WL,"Wavelength from Moon Plot"))
        self.log.debug(npArrayInfo(FL,"Flux from Moon Plot"))
        plt.semilogy(WL*1e6,FL,'m-',linestyle='-',label="Moon",zorder=0.5)
        
        plt.axis(axis)
        
        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("Flux (Photons)")
        plt.legend(loc=4)
        FileName = "%(Partials)s/Source-Spectrum%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
        
    
    @description("Plotting QE Spectrum")
    @depends("setup-sky")
    def plot_qe(self):
        """Plot sky spectrum"""
        WL = self.config["Instrument.wavelengths.values"]
        plt.figure()
        plt.clf()
        plt.title("QE Spectrum")
        WL,FL = self.qe.frame()(wavelengths=WL)
        plt.semilogy(WL*1e6,FL,'b.',linestyle='-',label="Quantum Efficiency")
        WL,FL = self.sky.frame("Atmosphere")(wavelengths=WL)
        plt.semilogy(WL*1e6,FL,'m.',linestyle='-',label="Extinction")
        
        ax = plt.gca()
        ax.yaxis.set_minor_formatter(
            LogFormatterTeXExponent(base=10,
             labelOnlyBase=False))
                     
        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("Flux (Fraction)")
        plt.legend(loc=8)
        FileName = "%(Partials)s/QE-Spectrum%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
    
    @description("Plotting Lenslet hexagons")
    @depends("setup-hexagons")
    def plot_hexagons(self):
        """Plot the lenslet Hexagons"""
        plt.figure()
        plt.clf()
        plt.title("Position of Lenslets")
        self.map_over_lenslets(lambda l:l.show_geometry(),color="cyan")
        plt.xlabel("x")
        plt.ylabel("y")
        FileName = "%(Partials)s/Lenslet-Hexagons%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
        
    
    @description("Plotting pixel squares")
    @depends("setup-source-pixels")
    def plot_pixels(self):
        """Plot the source pixels"""
        plt.figure()
        plt.clf()
        plt.title("Position of pixels")
        self.map_over_pixels(lambda p:p.show_geometry(),color="cyan")
        plt.xlabel("x")
        plt.ylabel("y")
        FileName = "%(Partials)s/Pixel-Geometry%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
    
    @description("Plotting invalid pixel squares")
    @depends("setup-source-pixels")
    def plot_invalid_pixels(self):
        """Plot invalid shapes"""
        plt.figure()
        plt.clf()
        plt.title("Position of invalid pixels")
        self.map_over_pixels(self._plot_invalid_pixels,color="cyan")
        plt.xlabel("x")
        plt.ylabel("y")
        FileName = "%(Partials)s/Pixel-Invalid-Geometry%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
    
    @ignore
    def _plot_invalid_pixels(self,pixel):
        """docstring for _plot_invalid_hexagon"""
        if not pixel.shape.is_valid:
            pixel.show_geometry()
    
    
    @description("Plotting invalid hexagons")
    @depends("setup-hexagons")
    def plot_invalid_hexagons(self):
        """Plot invalid shapes"""
        plt.figure()
        plt.clf()
        plt.title("Position of invalid hexagons")
        self.map_over_lenslets(self._plot_invalid_hexagon,color="cyan")
        plt.xlabel("x")
        plt.ylabel("y")
        FileName = "%(Partials)s/Lenslet-Invalid-Geometry%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
    
    @ignore
    def _plot_invalid_hexagon(self,lenslet):
        """docstring for _plot_invalid_hexagon"""
        if not lenslet.shape.is_valid:
            lenslet.show_geometry()
        
    
    @description("Plotting Lenslet-plane geometry")
    @depends("setup-source-pixels","setup-hexagons")
    def plot_geometry(self):
        """Plot all of the geomoetry stacked"""
        plt.figure()
        plt.clf()
        plt.title("Lenslet Plane Geometry")
        self.map_over_lenslets(lambda l:l.show_geometry(color="#cccc00"),color="cyan")
        self.map_over_pixels(lambda p:p.show_geometry(color="#cc00cc"),color="cyan")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.axes().set_aspect('equal')
        FileName = "%(Partials)s/System-Geometry%(fmt)s" % dict(fmt=self.config["Plots.format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
    
    
    @description("Plotting resample matrix")
    @depends("geometric-resample")
    def plot_resample(self):
        """Show the geometric resample"""
        self.map_over_pixels(self._show_resample,color="cyan")
    
    @ignore
    def _show_resample(self,pixel):
        """docstring for _show_resample"""
        plt.figure()
        plt.clf()
        plt.title("Resample for pixel %g" % pixel.num)
        pixel.show_geometry()
        self.map_over_lenslets(lambda l:self._show_lenslet_resample(l,pixel),color=False)
        # plt.colorbar()
        
        FileName = "%(Partials)s/System-Geometry-%(pixel)g%(fmt)s" % dict(pixel=pixel.num,fmt=self.config["Plots.format"],**self.config["Dirs"])
        
        plt.savefig(FileName)
        plt.clf()
    
    
    @description("Plotting lenslet ellipse sizes")
    @depends("setup-lenslets","dispersion")
    def plot_ellipses(self):
        """docstring for plot_ellipses"""
        plt.clf()
        self.map_over_lenslets(lambda l: l.plot_ellipses(),color="cyan")
        plt.title("Major Axis Size")
        plt.xlabel("Wavelength ($\mu m$)")
        plt.ylabel("$\Delta$-position ($px$)")
        plt.axis(expandLim(plt.axis()))
        plt.savefig("%(Partials)s/Lenslets-WL-dy%(ext)s" % dict(ext=self.config["Plots.format"],**self.config["Dirs"]))
        
    
    @description("Plotting lenslet ellipse sizes")
    @depends("setup-lenslets","dispersion")
    def plot_rotation(self):
        """Plot rotation"""
        plt.clf()
        self.map_over_lenslets(lambda l: l.plot_rotation(),color="cyan")
        plt.title("Major Axis Rotation")
        plt.xlabel("Wavelength ($\mu m$)")
        plt.ylabel("Rotation (Degrees)")
        plt.axis(expandLim(plt.axis()))
        plt.savefig("%(Partials)s/Lenslets-WL-dr%(ext)s" % dict(ext=self.config["Plots.format"],**self.config["Dirs"]))
        
    
    @description("Plotting PSF Kernels")
    @depends("setup-caches","setup-config")
    def plot_kernel_partials(self):
        """Plots the kernel data partials"""
        self.log.debug("Generating Kernel Plots and Images")
        major = self.config["Instrument.Tel.radius.px"] * self.config["Instrument.density"] * 1.2
        minor = self.config["Instrument.Tel.radius.px"] * self.config["Instrument.density"]
        ETEL = self.get_tel_kern(major,minor)
        ECONV = sp.signal.convolve(self.Caches["PSF"],ETEL,mode='same')
        plt.clf()
        plt.imshow(ETEL,interpolation='nearest')
        plt.title("Telescope Image (Ellipse)")
        plt.colorbar()
        plt.savefig("%s/Instrument-ETEL-Kernel%s" % (self.config["Dirs.Partials"],self.config["Plots.format"]))
        plt.clf()
        plt.clf()
        plt.imshow(ECONV,interpolation='nearest')
        plt.title("Convolved ETel + PSF Image")
        plt.colorbar()
        plt.savefig("%s/Instrument-ECONV-Kernel%s" % (self.config["Dirs.Partials"],self.config["Plots.format"]))
        plt.clf()
        
        
        plt.clf()
        plt.imshow(self.Caches["TEL"],interpolation='nearest')
        plt.title("Telescope Image")
        plt.colorbar()
        plt.savefig("%s/Instrument-TEL-Kernel%s" % (self.config["Dirs.Partials"],self.config["Plots.format"]))
        plt.clf()
        plt.imshow(self.Caches["PSF"],interpolation='nearest')
        plt.title("PSF Image")
        plt.colorbar()
        plt.savefig("%s/Instrument-PSF-Kernel%s" % (self.config["Dirs.Partials"],self.config["Plots.format"]))
        plt.clf()
        plt.imshow(self.Caches["CONV"],interpolation='nearest')
        plt.title("Convolved Tel + PSF Image")
        plt.colorbar()
        plt.savefig("%s/Instrument-FIN-Kernel%s" % (self.config["Dirs.Partials"],self.config["Plots.format"]))
        plt.clf()
        
        
    @ignore
    def _show_lenslet_resample(self,lenslet,pixel):
        """docstring for _show_lenslet_resample"""
        lenslet.show_geometry(color=pixel.get_color(lenslet.idx))
    
    
    @description("Writing resample Matrix")
    @depends("geometric-resample")
    def write_resample(self):
        """Write the resample matrix to file"""
        with open("%(Partials)s/Resample-info.dat" % dict(**self.config["Dirs"]),"w") as infostream:
            with open("%(Partials)s/Resample.dat" % dict(**self.config["Dirs"]),"w") as rawstream:
                infostream.write("# Resample Matrix \n")
                self.map_over_lenslets(lambda l: self._write_resample(l,infostream,rawstream),color="cyan")
    
    @ignore
    def _write_resample(self,lenslet,streama,streamb):
        """Write the resampleing matrix"""
        string = "%(lenslet)d %(info)s %(spec)s\n" % { 'lenslet': lenslet.num, 'info': npArrayInfo(lenslet.pixelValues), 'array': lenslet.pixelValues , 'spec' : str(lenslet.spectrum)}
        streama.write(string)
        np.savetxt(streamb,lenslet.pixelValues)
        
    
    #######################
    ## Mapping Functions ##
    #######################
    
    @ignore
    def map_over_lenslets(self,function,exceptions=True,color="green"):
        """Maps a given function to operate on each lenslet, and displays a progress bar along the way."""
        collection = self.lenslets.values()
        self.map_over_collection(function,lambda l:l.num,collection,exceptions,color)
        
    
    @ignore
    def map_over_pixels(self,function,exceptions=True,color="green"):
        """Maps some function over a bunch of source pixels"""
        collection = self.SourcePixels
        self.map_over_collection(function,lambda p:p.num,collection,exceptions,color)        
        
    
    ###################
    ## Image KERNELS ##
    ###################
    
    
    @ignore
    def generate_poisson_noise(self,label=None,lam=2.0,num=1):
        """Generates a poisson noise mask, saving to this object"""
        distribution = np.random.poisson
        shape = (self.config["Instrument.ccd.size.px"],self.config["Instrument.ccd.size.px"])
        if label == None:
            label = "Poisson Noise Mask (%2g)" % (lam)
        arguments = (lam,shape)
        noise = np.zeros(shape)
        for i in range(num):
            noise += distribution(*arguments)
        self.save(noise,label)
    
    @ignore
    def get_conv(self,wavelength,a=None,b=None,rot=None):
        """Return a PSF for a given wavelength in the system"""
        if a and b:
            ai = int(a)
            bi = int(b)
            ri = int(rot * 180.0 / np.pi)
            if ai not in self.ellipses:
                self.ellipses[ai] = {}
            if bi not in self.ellipses[ai]:
                self.ellipses[ai][bi] = {}
            if ri not in self.ellipses[ai][bi]:
                ETEL = self.get_tel_kern(a,b)
                self.CONV = sp.signal.fftconvolve(self.Caches["PSF"],ETEL,mode='same')
                self.ellipses[ai][bi][ri] = self.CONV
            else:
                self.CONV = self.ellipses[ai][bi][ri]
        if not hasattr(self,"found"):
            self.found = True
            self.CONV = self.Caches["CONV"] 
        return self.CONV
    
    @ignore
    def psf_kern(self,filename,size=0,truncate=False,header_lines=18):
        """Generates a PSF Kernel from a file with micron-encircled energy conversions. The file should have two columns, first, microns from the center of the PSF, and second, the fraction of encircled energy at that distance from the PSF.
        
        The calculation is then carried out using a spline fit to this data. From the spline fit, the function returns the first derivative of the encircled energy at each point. This in effect is the amount of energy at each point. These values are then normalized, to create a PSF mask for the instrument.
        
        The `size` parameter specifies the size of the kernel to use. If the size is greater than the encircled energy data, then a larger figure will be returned. If `size` is smaller than the encircled energy data, it will return an image the size of the encircled energy data, unless the `truncate` parameter is set to `true`.
        
        The `header_lines` parameter defaults to 18, which works with Zemax encircled energy output."""
        
        uM,FR = np.genfromtxt(filename,skip_header=header_lines).T
        # Convert microns to milimeters, then pixels, then dense pixels
        PX = uM * 1e-3 * self.config["Instrument.convert.mmtopx"] * self.config["Instrument.density"]
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
    
    @ignore
    def ellipse_kern(self,major,minor,alpha=0,size=0,sizey=False,normalize=False):
        """docstring for elipse_kern"""
        size /= 2
        sizey /= 2
        
        if size < sizey:
            size = sizey
        if size < major:
            size = int(major) + 1
        sizey = size
        
        major = float(major)
        minor = float(minor)
        alpha = float(alpha)
        
        x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
        
        d = np.sqrt(((x * np.cos(alpha) + y * np.sin(alpha))/minor)**2.0 + ((x * np.sin(alpha) + y * np.cos(alpha))/major)**2.0)
        
        v = (d <= 1).astype(np.float)
        if normalize:
            return v / np.sum(v)
        else:
            return v
        
        
    
    @ignore
    def circle_kern(self,radius,size=0,sizey=0,normalize=False):
        """Generate a Circle Kernel for modeling the \"Image of the Telescope\". The radius should be set in array units.
        
        `size` will determine the size of the array image, unless `size` is less than `radius`, in which case the image will be automatically increased to fit the entire circle.
        
        `normalize` controls whether the data is normalized or not. If it is not normalized, the data will have only 1.0 and 0.0 values, where 1.0 is within the radius, and 0.0 is outside the raidus."""
        size /= 2
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
    
    @ignore
    def gauss_kern(self,stdev,size=0,stdevy=None,sizey=None,enlarge=True,normalize=True):
        """ Returns a normalized 2D gaussian kernel array for convolutions.
        
        `stdev` is the standard deviation in the x-direction. If the `stdevy` keyword is not set, then it will be used as the standard deviation in the y-direction as well.
        
        `size` will determine the size of the returned image unless `size` is less than `stdev**2`.
        
        Results from this function are always normalized.
        """
        if size < (stdev*1.5) and enlarge:
            size = np.int(stdev*1.5)
        else:
            size = np.int(size)
        if not stdevy:
            stdevy = stdev
        if not sizey:
            sizey = size
        if sizey < (stdevy*1.5) and enlarge:
            sizey = np.int(stdevy*1.5)
        else:
            sizey = np.int(sizey)
        
        x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
        g = np.exp(-(x**2.0/np.float(2* (stdev**2.0))+y**2.0/np.float(2 * (stdevy**2.0))))
        
        if normalize:
            return g / g.sum()
        else:
            return g
    
    @ignore
    def get_tel_kern(self,major=None,minor=None):
        """Returns the telescope kernel. This kernel is built by creating a circle mask for the size of the telescope mirror, and then subtracting a telescope obscuration from the center of the mirror image. The values for all of these items are set in the configuration file."""
        if major or minor:
            
            
            
            TELIMG = self.ellipse_kern( major, minor )
            center = self.ellipse_kern( major * self.config["Instrument.Tel.obsc.ratio"], minor * self.config["Instrument.Tel.obsc.ratio"], *TELIMG.shape )
        else:
            TELIMG = self.circle_kern( self.config["Instrument.Tel.radius.px"] * self.config["Instrument.density"] )
            center = self.circle_kern( self.config["Instrument.Tel.obsc.px"] * self.config["Instrument.density"] ,
                *TELIMG.shape )
        TELIMG -= center
        TELIMG = TELIMG / np.sum(TELIMG)
        self.log.debug(npArrayInfo(TELIMG,"TELIMG"))
        return TELIMG
    
    @ignore
    def get_psf_kern(self):
        """Returns the PSF Kernel. The function first tries to read the encircled energy file. In this case, if a `psf_size` is set in the instrument configuration, this value will be used to truncate the size of the encircled energy function. If the encircled energy function cannot be loaded, the system will fall back on to a gaussian psf as configured by the instrument."""
        if self.config["Instrument.PSF.size.px"] != 0:
            size = self.config["Instrument.PSF.size.px"] * self.config["Instrument.density"]
            truncate = True
        else:
            size = 0
            truncate = False
        try:
            PSFIMG = self.psf_kern(self.dir_filename("Data",self.config["Instrument.files.encircledenergy"]),size,truncate)
        except IOError as e:
            self.log.warning("Could not access encircled energy file: %s" % e)
            PSFIMG = self.gauss_kern( (self.config["Instrument.PSF.stdev.px"] * self.config["Instrument.density"]) )
        else:
            self.log.debug("Loaded Encircled Energy from %s" % self.config["Instrument.files.encircledenergy"])
        return PSFIMG
    
    @ignore
    @description("Setting up configuration values")
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
        self.config.load()
        self.astrologger.configure(configFile = self.config["Configurations.This"])
        self.astrologger.start()
        self.astrologger.useConsole(False)
        
        if "calc" not in self.config["Instrument.convert"]:
            if "mmtopx" not in self.config["Instrument.convert"] and "pxtomm" in self.config["Instrument.convert"]:
                self.config["Instrument.convert.mmtopx"] = 1.0 / self.config["Instrument.convert.pxtomm"]
                self.config["Instrument.convert.calc"] = True
            else:
                self.config["Instrument.convert.pxtomm"] = 1.0 / self.config["Instrument.convert.mmtopx"]
                self.config["Instrument.convert.calc"] = True
        
        self.config["Instrument"] = self._setUnits(self.config["Instrument"],None)
        
        self.config["Instrument.image.size.px"] = np.round( self.config["Instrument.image.size.px"] , 0 )
                
        wl,r = self.get_resolution_spectrum(self.config["Instrument.wavelengths.min"],self.config["Instrument.wavelengths.max"],self.config["Instrument.wavelengths.resolution"])
        
        self.config["Instrument.wavelengths.values"] = wl
        self.config["Instrument.wavelengths.resolutions"] = r
        
        sys.setrecursionlimit(10000000)
    
    @ignore
    def get_resolution_spectrum(self,minwl,maxwl,resolution):
        """docstring for get_resolution_spectrum"""
        
        dwl = [minwl]
        new_wl = minwl
        while new_wl <= maxwl:
            new_wl += (new_wl / resolution)
            dwl += [new_wl]
            
        dense_wavelengths = np.array(dwl)
        dense_resolution = dense_wavelengths[:-1] / np.diff(dense_wavelengths)
        dense_wavelengths = dense_wavelengths[:-1]
        return dense_wavelengths, dense_resolution
    
    @ignore
    def _setUnits(self,r,parent):
        """docstring for _setUnits"""
        for k in r.keys():
            v = r[k]
            if isinstance(v, collections.Mapping):
                r[k] = self._setUnits(r[k],k)
            elif k == "mm":
                if ("calc" in r) and ("px" in r):
                    pass
                elif ("calc" not in r):
                    r["px"] = v * self.config["Instrument.convert.mmtopx"]
                    r["calc"] = True
                else:
                    self.log.warning("Value for %s set in both px and mm." % parent)
            elif k == "px":
                if ("calc" in r) and ("mm" in r):
                    pass
                elif ("calc" not in r):
                    r["mm"] = v * self.config["Instrument.convert.pxtomm"]
                    r["calc"] = True
                else:
                    self.log.warning("Value for %s set in both px and mm." % parent)
        return r

def run():
    SIM = SEDSimulator()
    SIM.run()
    
if __name__ == '__main__':    
    run()
