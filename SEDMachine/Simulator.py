#!/usr/bin/env python
# 
#  Simulator.py
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2012-02-08.
#  Copyright 2012 Alexander Rudy. All rights reserved.
#  Version 0.3.0
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
import gc

import AstroObject
from AstroObject.AstroSimulator import Simulator
from AstroObject.AstroCache import *
from AstroObject.AstroConfig import *
from AstroObject.AstroSpectra import SpectraObject
from AstroObject.AstroImage import ImageObject,ImageFrame
from AstroObject.AnalyticSpectra import BlackBodySpectrum, AnalyticSpectrum, FlatSpectrum, ResampledSpectrum, FLambdaSpectrum, InterpolatedSpectrum
from AstroObject.Utilities import *

from Lenslet import *


class SEDSimulator(Simulator,ImageObject):
    """A simulator for the SED Machine"""
    def __init__(self):
        super(SEDSimulator, self).__init__(name="SEDMachine")
        self.debug = False
        self.dataClasses = [SubImage]
        self.lenslets = []
        self.config.merge(self.basics)
        self.config.merge({"Instrument":self.instrument,"Caches":self.caches,"Source":self.source,"Observation":self.observation})
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
            "FrameLabel": "Final",
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
                'level': logging.INFO, 
                'enable': True, 
                'format': '%(levelname)-8s... %(message)s'
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
    
    observation = {
        'exposure' : 1200,
        'number' : 3,
        'airmass' : 1,
        
    }
    
    
    instrument = {
        'files': {
            'dispersion': 'Data/dispersion_12-10-2011.txt',
            'encircledenergy': 'Data/encircled_energy_4nov11.TXT',
            'lenslets': 'Data/xy_17nov2011_v57.TXT',
        },
        'camera' : "PI",
        'convert': {
            'pxtomm': 0.0135 }, 
        'density': 5, 
        'tel_obsc': {'px': 0.2}, 
        'plot': False, 
        'ccd_size': {'px': 2048}, 
        'padding': 5, 
        'psf_stdev': {'px': 1.0}, 
        'bias': 20, 
        'psf_size': { 'px': 2.4}, 
        'image_size': { 'mm': 40.0}, 
        'tel_radii': { 'px': 1.2},
        'tel_area' : 18242. * 0.9,
        'gain': 5,
        'eADU' : 3.802,
        'lenslets' : {
              'radius' : 0.25e-2,
              'rotation' : 20.0,
        },
        'Sky' : {
            'Use' : "TurnroseSKY",
            'Files' : {
                'Massey' : "SEDSpec2/MasseySky.fits",
                'Quimby' : "SEDSpec2/Quimby.fits",
                'HansuchikUVES' : "SEDSpec2/HansuchikUVES.fits",
                'TurnroseSKY' : "SEDSpec2/TurnroseSKY.fits",
            },
        },
        'Thpt' : {
            'File' : "SEDSpec2/Data/thpt.npy",
            'Type' : "prism_pi",
        },
        'Moon' : {
            'Phase' : 1,
        },
        
    }
    
    source = {
        'Filename' : "Data/SNIa.R1000.dat",
        'CubeName' : "Data/CUBE.fits",
        'PreAmp' : 100.0,
        'PXSize' : { 'mm' : 0.5 },
        'Rotation' : 0,
    }
    
    def setup_stages(self):
        """Sets up all simulator stages"""
        self.registerConfigOpts("D",{"Lenslets":{"start":2100,"number":50},"Debug":True},help="Debug, Limit lenslets (50,start from 2100)")
        self.registerConfigOpts("S",{"Lenslets":{"start":2100,"number":5},"Debug":True},help="Debug, Limit lenslets (5,start from 2100)")
        self.registerConfigOpts("T",{"Lenslets":{"start":2100,"number":50},"Debug":False},help="Limit lenslets (50,start from 2100)")
        self.registerConfigOpts("M",{"Lenslets":{"start":1000,"number":500},"Debug":True},help="Limit lenslets (500,start from 1000)")
        self.registerConfigOpts("A",{"Lenslets":{"start":2100,"number":1},"Debug":True},help="Debug, Single lenslets (start from 2100)")
        
        
        
        self.registerStage(self.setup_caches,"setup-caches",help=False,description="Setting up caches")
        self.registerStage(self.setup_configuration,"setup-config",help=False,description="Setting up dynamic configuration")
        self.registerStage(self.setup_constants,"setup-constants",help=False,description="Setting up physical constants")
        self.registerStage(self.setup_cameras,"setup-cameras",help=False,description="Setting up Cameras")
        self.registerStage(self.setup_lenslets,"setup-lenslets",help=False,description="Setting up lenslets",dependencies=["setup-config"])
        self.registerStage(self.setup_hexagons,"setup-hexagons",help=False,description="Setting up lenslet hexagons",dependencies=["setup-lenslets"])
        self.registerStage(self.setup_blank,"setup-blank",help=False,description="Creating blank image",dependencies=["setup-config"])
        self.registerStage(self.setup_source,"setup-source",help=False,description="Creating source spectrum objects",dependencies=["setup-config","setup-constants"])
        self.registerStage(self.setup_simple_source,"setup-source-simple",help=False,description="Creating simple source spectrum object",dependencies=["setup-config","setup-constants"],include=False,replaces=["setup-source"])
        self.registerStage(self.setup_noise,"setup-noise",help=False,description="Setting up Dark/Bias frames",dependencies=["setup-config","setup-cameras"])
        self.registerStage(self.setup_sky,"setup-sky",help=False,description="Setting up Sky spectrum object",dependencies=["setup-config","setup-constants"])
        self.registerStage(self.flat_source,"flat-source",help="Make a constant value source",description="Replacing default source with a flat one.",include=False,replaces=["setup-source"])
        self.registerStage(None,"simple-source",help="Use a simple, centered source object",description="Replacing default source with a simple one",include=False,dependencies=["setup-source-simple"])
        self.registerStage(None,"setup",help="System Setup",description="Set up simulator",
            dependencies=["setup-caches","setup-lenslets","setup-hexagons","setup-blank","setup-source","setup-noise","setup-constants","setup-sky","setup-cameras"],
            )
        
        self.registerStage(self.apply_qe,"apply-qe",help=False,description="Applying Quantum Efficiency Functions",dependencies=["setup-source","setup-sky"])
        
        self.registerStage(self.use_sky,"use-sky",help="Use only sky spectrum",description="Using only Sky spectrum",dependencies=["setup-sky","apply-qe"],include=False)
        
        self.registerStage(self.plot_source,"plot-source",help="Plot sky spectrum",description="Plotting Source Spectrum",include=False,dependencies=["setup-sky","setup-source","apply-qe"])
        self.registerStage(self.plot_sky,"plot-sky",help="Plot sky spectrum",description="Plotting Sky Spectrum",include=False,dependencies=["setup-sky","apply-qe"])
        self.registerStage(self.plot_qe,"plot-qe",help="Plot QE spectrum",description="Plotting QE Spectrum",include=False,dependencies=["setup-sky"])
        self.registerStage(self.plot_lenslet_data,"plot-lenslet-xy",help="Plot Lenslets",description="Plotting lenslet positions",include=False,dependencies=["setup-lenslets"])
        
        self.registerStage(self.lenslet_dispersion,"dispersion",help="Calculate dispersion",description="Calculating dispersion for each lenslet",dependencies=["setup-lenslets","setup-caches"])
        self.registerStage(self.lenslet_trace,"trace",help="Trace Lenslets",description="Tracing lenslet dispersion",dependencies=["dispersion","setup-caches","apply-qe"])
        self.registerStage(self.lenslet_place,"place",help="Place Subimages",description="Placing lenslet spectra",dependencies=["trace"])
        
        self.registerStage(self.plot_dispersion_data,"plot-dispersion",help=False,description="Plotting dispersion for each lenslet",dependencies=["dispersion"],include=False)
        self.registerStage(self.plot_trace_data,"plot-trace",help=False,description="Plotting trace data for each lenslet",dependencies=["trace"],include=False)
        self.registerStage(self.plot_spectrum_data,"plot-spectrum",help=False,description="Plotting spectral data for each lenslet",dependencies=["trace"],include=False)
        self.registerStage(None,"plot-lenslets",help="Plot data about each lenslet",description="Plotting data about each lenslet",include=False,dependencies=["plot-dispersion","plot-trace","plot-spectrum"])
        
        self.registerStage(self.image_merge,"merge-cached",help="Merge Cached Subimages",description="Merging subimages",dependencies=["setup-blank","setup-lenslets"])
        self.registerStage(None,"merge",help="Merge Subimages",description="Merging subimage into master image",dependencies=["place","merge-cached"])
        
        self.registerStage(self.ccd_crop,"crop",help="Crop Final Image",description="Cropping image to CCD size",dependencies=["setup-blank"])
        self.registerStage(self.apply_noise,"add-noise",help="Add Dark/Bias noise to image",description="Adding Dark/Bias noise",dependencies=["crop","setup-noise"])
        self.registerStage(self.save_file,"save",help="Save image to file",description="Saving image to disk",dependencies=["setup-blank"])
        self.registerStage(None,"cached-only",help="Use cached subimages to construct final image",description="Building image from caches",dependencies=["merge-cached","crop","add-noise","save"],include=False)
        
        self.registerStage(None,"plot",help="Create all plots",description="Plotting everything",dependencies=["plot-lenslet-xy","plot-lenslets","plot-sky","plot-qe","plot-source"],include=False)
        
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
    
    def setup_constants(self):
        """Establish Physical Constants"""
        self.const = StructuredConfiguration()
        self.const.setFile("const","SED.const.config.yaml")
        self.const["hc"] = 1.98644521e-8 # erg angstrom
        
        
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
       FileName = "%(Partials)s/%(name)s%(ext)s" % dict(name="Lenslets-raw",ext=".dat",**self.config["Dirs"])
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
       
       
       # Central Lenslet Output
       FileName = "%(Partials)s/%(name)s%(ext)s" % dict(name="center-raw",ext=".dat",**self.config["Dirs"])
       with open(FileName,'w') as stream:
           lenslet = self.lenslets[ix[cntix]]
           stream.write(lenslet.introspect())
       
       
       if "start" in self.config["Lenslets"]:
           self.lensletIndex = self.lensletIndex[self.config["Lenslets"]["start"]:]
       if "number" in self.config["Lenslets"]:
           self.lensletIndex = self.lensletIndex[:self.config["Lenslets"]["number"]]
       self.total = len(self.lensletIndex)
       self.lenslets = {x:self.lenslets[x] for x in self.lensletIndex}
    
    def setup_blank(self):
        """Establish a blank Image"""
        self.save(np.zeros((self.config["Instrument"]["image_size"]["px"],self.config["Instrument"]["image_size"]["px"])).astype(np.int16),"Blank")

    def setup_source(self):
        """Sets up a uniform source file based spectrum"""
        
        
        Source = ImageObject()
        Source.read(self.config["Source"]["CubeName"])
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
                self.SourcePixels.append(SourcePixel(i,j,data=data[i,j],label="Source Pixel %d,%d" % (i,j)))
                progress = int((finished/float(total)) * 100)
                finished += 1
                PBar.render(progress,"L:%4d,%4d %6d/%-6d" % (i,j,finished,total))
        
        PBar.render(1,"L:%4s,%4s %6d/%-6d" % ("Done","",finished,total))
        
        
    def setup_simple_source(self):
        """docstring for setup_simple_source"""        
        WL,FL = np.genfromtxt(self.config["Source"]["Filename"]).T
        FL /= self.const["hc"] / WL
        FL *= 1e10 #Spectrum was per Angstrom, should now be per Meter
        WL *= 1e-10
        self.Spectrum = FLambdaSpectrum(np.array([WL,FL]),self.config["Source"]["Filename"])
        self.Original = FLambdaSpectrum(np.array([WL,FL]),self.config["Source"]["Filename"])
        self.SourcePixels = [SourcePixel(0,0,data=np.array([WL,FL]),label="Source Pixel",config=self.config)]
    
    
    def setup_cameras(self):
        """Set up camera configuration values"""
        # Camera Data (from sim_pdr.py by Nick)
        self.cameras = {
        	"PI": {"DQEs": np.array([
        		(2000, 0),	
        		(3000, 0.01),
        		(3500, .20),
        		(4000, .60),
        		(4500, .82),
        		(5000, .90),
        		(5500, .93),
        		(6000, .93),
        		(7000, .93),
        		(7500, .88),
        		(8000, .73),
        		(8500, .55),
        		(9000, .33),
        		(10000, .08),
        		(10500, 0.02),
        		(11000, 0)
        	]),
        	"RN" : 5.,
        	"DC":  0.006,
        	"readtime": 37},
        	"Andor": { # Midband
        		"DQEs": np.array([
        		(2500, .05),
        		(3500, .18),
        		(4500, .75),
        		(5000, .9),
        		(5500, .92),
        		(6500, .91),
        		(7500, .79),
        		(8500, .48),
        		(9500, .13),
        		(10500, .02),
        		(11000, 0)
        	]), 
        	"RN": 4,
        	"DC": 0.0004,
        	"readtime":  82},
        	"E2V": {"DQEs" : np.array([
        		(3000, .1),
        		(3500, .3),
        		(4500, .8),
        		(5500, .8),
        		(6500, .78),
        		(7500, .7),
        		(8500, .4),
        		(9500, .13),
        		(10500, .02),
        		(11000, 0)]),
        	"RN": 3.3,
        	"DC": 0.006,
        	"readtime": 37},
        }


        self.cameras["PI-fast"] = self.cameras["PI"]
        self.cameras["PI-fast"]["RN"] = 12
        self.cameras["PI-fast"]["readtime"] = 2.265

        self.cameras["Andor-fast"] = self.cameras["Andor"]
        self.cameras["Andor-fast"]["RN"] = 11.7
        self.cameras["Andor-fast"]["readtime"] = 1.398
        
    
    def setup_sky(self):
        """Setup sky spectrum information"""

        
        # Sky Data (From sim_pdr.py by Nick, regenerated using SEDSpec2 module's make_files.py script)
        # Each sky spectrum is saved in a FITS file for easy recall as a spectrum object.
        self.SKYData = SpectraObject()
        for filename in self.config["Instrument"]["Sky"]["Files"].values():
            self.SKYData.read(filename)
        
        # Moon phase adjustments. These moon phase attenuation values are for different filter bands.
        # The intermediate wavelengths are accounted for using a polyfit
        # the result is a function which takes phase and wavelength, and outputs an attenuation...
        
        # See derivation on pg 83 of SED NB 1 (20 July 2011)
        moon_phase = np.array([0., 0.08, 0.16, 0.24, 0.32, 0.40, 0.50])
        moon_g = np.array([2e-17, 2.1e-17, 2.15e-17, 2.3e-17, 5.3e-17, 1.7e-16, 3.2e-16])
        moon_r = np.array([2.3e-17,2.3e-17,2.3e-17,3.3e-17,3.5e-17,8.3e-17,1.3e-16])
        moon_i = np.array([2.8e-17,3.0e-17,3.0e-17,3.3e-17,3.8e-17,7.0e-17,9.0e-17])

        sky_ls = np.array([4868., 6290., 7706., 10000])

        self.moon_specs = []
        for i in xrange(len(moon_phase)):
            gm = moon_g[i]-moon_g[0]
            rm = moon_r[i]-moon_r[0]
            im = moon_i[i]-moon_i[0]
            zm = im
            fluxes = np.array([gm, rm, im, zm])
            # fluxes /= self.const["hc"] / sky_ls
            moon_spec = InterpolatedSpectrum(np.array([sky_ls,fluxes]),"Moon Phase %s" % i)
            moon_spec.func = np.poly1d(np.polyfit(sky_ls,fluxes, 2))
            self.moon_specs.append(moon_spec)
        
        # Throughputs are generated from Nick's simulation scripts in throughput.py
        # They are simply re-read here.
        thpts = np.load(self.config["Instrument"]["Thpt"]["File"])[0]
        self.qe = {}
        WL = thpts["lambda"]* 1e-10
        self.qe["prism_pi"] = InterpolatedSpectrum(np.array([WL, thpts["thpt-prism-PI"]]),"PI Prism")
        self.qe["prism_andor"] = InterpolatedSpectrum(np.array([WL, thpts["thpt-prism-Andor"]]),"Andor Prism")
        self.qe["grating"] = InterpolatedSpectrum(np.array([WL, thpts["thpt-grating"]]),"Grating")
        
        # This calculation fixes the units of the TurnroseSKY values
        # I'm not sure what these units are doing, but we will leave them here for now.
        WL,FL = self.SKYData.data(self.config["Instrument"]["Sky"]["Use"])
        FL *= 1e-18 * 3 # NICK! I NEED THIS EXPLAINED!
        FL += self.moon_specs[self.config["Instrument"]["Moon"]["Phase"]](wavelengths=WL)[1]
        FL /= self.const["hc"] / WL
        FL *= 1e10 #Spectrum was per Angstrom, should now be per Meter
        WL *= 1e-10
        
        self.Sky_Original = InterpolatedSpectrum(np.array([WL,FL]),"SkySpectrum")
        self.SkySpectrum = FLambdaSpectrum(np.array([WL,FL]),"SkySpectrum")
        
    def setup_hexagons(self):
        """Make the lenslet hexagons"""
        self.map_over_lenslets(lambda l: l.make_hexagon(),color="green")
        
    def apply_qe(self):
        """Apply the instrument quantum efficiency"""
        self.SkySpectrum *= self.qe[self.config["Instrument"]["Thpt"]["Type"]] * self.config["Instrument"]["tel_area"] * self.config["Observation"]["exposure"]
        self.Spectrum *= self.qe[self.config["Instrument"]["Thpt"]["Type"]] * self.config["Instrument"]["tel_area"] * self.config["Observation"]["exposure"]
        self.Original *= self.qe[self.config["Instrument"]["Thpt"]["Type"]] * self.config["Instrument"]["tel_area"] * self.config["Observation"]["exposure"]
        self.Spectrum += self.SkySpectrum
        
        
    def use_sky(self):
        """Use the sky spectrum only"""
        self.Spectrum = self.SkySpectrum
        
    def flat_source(self):
        """Replace the default file-source with a flat spectrum"""
        self.Spectrum = FlatSpectrum(self.config["Source"]["Value"]) * self.config["Source"]["PreAmp"]

    def lenslet_dispersion(self):
        """Calculate the dispersion for each lenslet"""
        self.map_over_lenslets(lambda l: l.find_dispersion(),color="blue")
        
    def lenslet_trace(self):
        """Trace out each lenslet"""
        self.map_over_lenslets(lambda l: l.get_trace(self.Spectrum),color="blue")

    def lenslet_place(self):
        """Place each spectrum into the subimage"""
        self.map_over_lenslets(self._lenslet_place,color="yellow")
        
    def _lenslet_place(self,l):
        """docstring for _lenslet_place"""
        l.place_trace(self.get_psf)
        l.write_subimage()
        with open("%(Partials)s/LensletAudit.dat" % self.config["Dirs"],"a") as s:
            s.write("%s\n" % vars(l) )
        # gc.collect()
    
    def image_merge(self):
        """Merge subimages into master image"""
        self.save(self.frame("Blank"),self.config["Output"]["FrameLabel"])
        self.map_over_lenslets(self._lenslet_merge,color="yellow")
        
    def _lenslet_merge(self,lenslet):
        """Merge a single lenslet into the master image"""
        lenslet.read_subimage()
        lenslet.bin_subimage()
        self.place(lenslet.data(),lenslet.subcorner,self.config["Output"]["FrameLabel"])
        lenslet.clear(delete=True)
    
    def ccd_crop(self):
        """Crops the image to the appropriate ccd size"""
        x,y = self.center
        size = self.config["Instrument"]["ccd_size"]["px"] / 2.0
        self.crop(x,y,size,label=self.statename)
        
    def setup_noise(self):
        """Makes noise masks"""
        
        read_noise = np.sqrt(self.cameras[self.config["Instrument"]["camera"]]["RN"]**2 * self.config["Instrument"]["psf_size"]["px"] * self.config["Observation"]["number"] / self.config["Instrument"]["gain"])
        dark_noise = self.cameras[self.config["Instrument"]["camera"]]["DC"] * self.config["Observation"]["exposure"]
        
        self.generate_poisson_noise("Read",read_noise)
        
        self.generate_poisson_noise("Dark",dark_noise)
    
        
    def apply_noise(self):
        """Apply the noise masks to the target image label"""
        
        dark = self.data("Dark")
        bias = self.data("Read")
        
        data = self.data(self.config["Output"]["FrameLabel"])
        
        data += dark + bias
        
        self.remove(self.config["Output"]["FrameLabel"])
        self.save(data,self.config["Output"]["FrameLabel"])
    
    def save_file(self):
        """Saves the file"""
        self.Filename = "%(Images)s/%(label)s-%(date)s.%(fmt)s" % dict(label=self.config["Output"]["Label"],date=time.strftime("%Y-%m-%d"), fmt=self.config["Output"]["Format"], **self.config["Dirs"] )
        self.write(self.Filename,states=[self.statename],clobber=True)
        self.log.info("Wrote %s" % self.Filename)
    
    
    ################################
    ## IMAGE MANAGEMENT FUNCTIONS ##
    ################################
        
    def place(self,img,corner,label):
        """Place the given AstroObject.AnalyticSpectra.AnalyticSpectrum onto the SEDMachine Image"""
        
        xstart = corner[0]
        xend = xstart + img.shape[0]
        ystart = corner[1]
        yend = ystart + img.shape[1]
        data = self.data(label)
        
        if data.shape[0] < xend or data.shape[1] < yend:
            raise SEDLimits
        
        if xstart < 0 or ystart < 0:
            raise SEDLimits
        
        if xend < 0 or yend < 0:
            raise SEDLimits
        
        data[xstart:xend,ystart:yend] += img
        self.remove(label)
        self.save(data,label)    
        
    
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
    
    
    #######################
    ## DEBUGGING METHODS ##
    #######################
    
    def plot_dispersion_data(self):
        """Outputs dispersion debugging data"""
        self.map_over_lenslets(lambda l: l.plot_dispersion(),color="cyan")
    
    def plot_trace_data(self):
        """Outputs plots about each lenslet trace"""
        self.map_over_lenslets(lambda l: l.plot_trace(),color="cyan")
        
    def plot_spectrum_data(self):
        """Outputs plots about each lenslet trace"""
        self.map_over_lenslets(lambda l: l.plot_spectrum(),color="cyan")
    
    def plot_lenslet_data(self):
        """Outputs the lenslet data"""
        plt.figure()
        plt.clf()
        self.log.info("Plotting lenslet arc positions in CCD (x,y) space")
        FileName = "%(Partials)s/Lenslet-xy%(fmt)s" % dict(fmt=self.config["plot_format"],**self.config["Dirs"])
        self.map_over_lenslets(lambda l: plt.plot(l.xs,l.ys,linestyle='-'),color="cyan")
        plt.title("Lenslet x-y positions")
        plt.savefig(FileName)
        
        plt.clf()
        self.log.info("Plotting lenslet physical positions in mm space")
        FileName = "%(Partials)s/Lenslet-pxy%(fmt)s" % dict(fmt=self.config["plot_format"],**self.config["Dirs"])
        self.map_over_lenslets(lambda l: plt.plot(l.ps.T[0],l.ps.T[1],marker='.'),color="cyan")
            
            
        plt.title("Lenslet p-xy positions")
        plt.savefig(FileName)
        
        plt.clf()
    
    def plot_sky(self):
        """Plot sky spectrum"""
        wavelengths, flux = self.SKYData.data(self.config["Instrument"]["Sky"]["Use"])
        WL = wavelengths*1e-10
        DWL = np.diff(WL) 
        WL = WL[:-1]
        RS = np.ones(WL.shape) * 100
        plt.figure()
        plt.title("Resolution")
        plt.plot(WL*1e6,RS,'g.')
        FileName = "%(Partials)s/Sky-Spectrum-Res%(fmt)s" % dict(fmt=self.config["plot_format"],**self.config["Dirs"])
        plt.savefig(FileName)

        plt.clf()
        plt.title("Sky Spectrum")
        self.log.debug(npArrayInfo(WL,"Wavelength for Sky Plot"))
        WL,FL = self.SkySpectrum(wavelengths=WL,resolution=RS)
        self.log.debug(npArrayInfo(WL,"Wavelength from Sky Plot"))
        self.log.debug(npArrayInfo(FL,"Flux from Sky Plot"))
        plt.semilogy(WL*1e6,FL,'b.',linestyle='-')
        WL,FL = self.Sky_Original(wavelengths=WL,resolution=RS)
        self.log.debug(npArrayInfo(WL,"Wavelength from Sky-O Plot"))
        self.log.debug(npArrayInfo(FL,"Flux from Sky-O Plot"))
        plt.semilogy(WL*1e6,FL,'g.',linestyle='-')
        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("Flux (Photons)")
        FileName = "%(Partials)s/Sky-Spectrum%(fmt)s" % dict(fmt=self.config["plot_format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()

    def plot_source(self):
        """Plot the source spectrum"""
        WL = np.arange(3300,10000,100) * 1e-10
        DWL = np.diff(WL) 
        WL = WL[:-1]
        RS = np.ones(WL.shape) * 100
        plt.figure()
        plt.title("Resolution")
        plt.plot(WL*1e6,RS,'g.')
        FileName = "%(Partials)s/Source-Spectrum-Res%(fmt)s" % dict(fmt=self.config["plot_format"],**self.config["Dirs"])
        plt.savefig(FileName)

        plt.clf()
        plt.title("Source Spectrum")
        self.log.debug(npArrayInfo(WL,"Wavelength for Sky Plot"))
        WL,FL = self.Spectrum(wavelengths=WL,resolution=RS)
        self.log.debug(npArrayInfo(WL,"Wavelength from Sky Plot"))
        self.log.debug(npArrayInfo(FL,"Flux from Sky Plot"))
        plt.semilogy(WL*1e6,FL,'b.',linestyle='-')
        WL,FL = self.SkySpectrum(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'g.',linestyle='-')
        WL,FL = self.Original(wavelengths=WL,resolution=RS)
        plt.semilogy(WL*1e6,FL,'r.',linestyle='-')
        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("Flux (Photons)")
        FileName = "%(Partials)s/Source-Spectrum%(fmt)s" % dict(fmt=self.config["plot_format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
        

    def plot_qe(self):
        """Plot sky spectrum"""
        wavelengths, flux = self.SKYData.data(self.config["Instrument"]["Sky"]["Use"])
        WL = wavelengths[1:]*1e-10
        plt.figure()
        plt.clf()
        plt.title("QE Spectrum")
        WL,FL = self.qe[self.config["Instrument"]["Thpt"]["Type"]](wavelengths=WL)
        plt.semilogy(WL*1e6,FL,'bo')
        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("Flux (Fraction)")
        FileName = "%(Partials)s/QE-Spectrum%(fmt)s" % dict(fmt=self.config["plot_format"],**self.config["Dirs"])
        plt.savefig(FileName)
        plt.clf()
    
    
    def map_over_lenslets(self,function,exceptions=True,color="green"):
        """Maps a given function to operate on each lenslet, and displays a progress bar along the way."""
        if exceptions == True:
            exceptions = Exception
        collection = self.lenslets.values()
        self.progress = 0.0
        self.errors = 0
        self.bar = arpytools.progressbar.ProgressBar(color=color)
        self.log.useConsole(False)
        self.bar.render(0,"L:%4s %4d/%-4d" % ("",self.progress,self.total))
        map(lambda l:self._collection_map(l,function,exceptions,l.num),collection)
        self.bar.render(100,"L:%4s %4d/%-4d" % ("Done",self.progress,self.total))
        self.log.useConsole(True)
        if self.progress != self.total:
            self.log.warning("Progress and Total are different at end of loop: %d != %d" % (self.progress,self.total))
        if self.errors != 0:
            self.log.warning("Trapped %d errors" % self.errors)
        
    def map_over_pixels(self,function,exceptions=True,color="green"):
        """Maps some function over a bunch of source pixels"""
        if exceptions == True:
            exceptions = Exception
        collection = self.SourcePixels
        self.progress = 0.0
        self.errors = 0
        self.bar = arpytools.progressbar.ProgressBar(color=color)
        self.log.useConsole(False)
        self.bar.render(0,"L:%4s %4d/%-4d" % ("",self.progress,self.total))
        map(lambda l:self._collection_map(l,function,exceptions,l.num),collection)
        self.bar.render(100,"L:%4s %4d/%-4d" % ("Done",self.progress,self.total))
        self.log.useConsole(True)
        if self.progress != self.total:
            self.log.warning("Progress and Total are different at end of loop: %d != %d" % (self.progress,self.total))
        if self.errors != 0:
            self.log.warning("Trapped %d errors" % self.errors)
        
            
    def _collection_map(self,lenslet,function,exceptions,identity):
        """Maps something over a bunch of lenslets"""
        self.bar.render(int(self.progress/self.total * 100),"L:%4d %4d/%-4d" % (identity,self.progress,self.total))
        try:
            function(lenslet)
        except exceptions as e:
            self.log.error("Caught %s in %d" % (e.__class__.__name__,identity))
            self.log.error(str(e))
            self.errors += 1
            if self.config["Debug"]:
                raise
        finally:
            self.progress += 1.0
            self.bar.render(int(self.progress/self.total * 100),"L:%4d %4d/%-4d" % (identity,self.progress,self.total))
        
        
        
    ###################
    ## Image KERNELS ##
    ###################
    
    
    
    def generate_poisson_noise(self,label=None,lam=2.0):
        """Generates a poisson noise mask, saving to this object"""
        distribution = np.random.poisson
        shape = (self.config["Instrument"]["ccd_size"]["px"],self.config["Instrument"]["ccd_size"]["px"])
        if label == None:
            label = "Poisson Noise Mask (%2g)" % (lam)
        arguments = (lam,shape)
        noise = distribution(*arguments)
        self.save(noise,label)
    
    def get_psf(self,wavelength):
        """Return a PSF for a given wavelength in the system"""
        if not hasattr(self,"found"):
            self.found = True
            self.CONV = self.Caches.get("CONV") 
        return self.CONV
    
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

def run():
    SIM = SEDSimulator()
    SIM.run()
    
if __name__ == '__main__':    
    run()
