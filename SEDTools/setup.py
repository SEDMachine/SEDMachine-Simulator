#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  setup.py
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2012-03-29.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

import shutil,os

import numpy as np

from pkg_resources import resource_filename

from AstroObject.AstroSimulator import *
from AstroObject.AstroSpectra import SpectraObject
import AstroObject.Utilities as AOU

# Background Functions
def abmag_to_flambda(AB , lam):
    # Ab magnitude
    # Wavelength in angstrom
    c = 2.9979e18 # angstrom/s
    # return erg/s/cm^2/ang
    return 10**(-(AB + 48.6)/2.5)*c/lam**2


class SetupAgent(Simulator):
    """An agent for setting up directories for simulations."""
    def __init__(self):
        super(SetupAgent, self).__init__(commandLine=True)
        self.config.load(resource_filename("SEDMachine","SED.main.config.default.yaml"))
        self.config.setFile("Main")
        self.config.load()
        self.config.update({"Default":"*all"})
        self.collect()
        self.registerStage(None,"make-specs",dependencies=["make-massey","make-hansuchik","make-turnrose","make-quimby","make-atmosphere","make-palsky"],include=True,help="Setup spectral FITS files")
     
    @include()
    @description("Setting up directories")
    @help("Make required Directories")   
    def setup_dirs(self):
        """Ensure that the required directories exist."""
        new_dirs = []
        for key,directory in self.config["Dirs"].iteritems():
            try:
                os.mkdir(directory)
            except os.error, e:
                self.log.debug("Directory %s already exists." % directory)
            else:
                self.log.debug("Created directory %s" % directory)
                new_dirs.append(directory)
        if len(new_dirs) > 0:
            self.log.info("Created Directories %r" % new_dirs)
        else:
            self.log.info("No Directories needed to be created.")
    
    @description("Collecting static (.dat) files")   
    @help("Populate data directory with included static data")
    @depends("setup-dirs")
    def get_data(self):
        """Copy static data files to the data directory."""
        dataDir = self.config["Dirs"]["Data"] + "/"
        sourceDir = resource_filename(__name__,"Data/")
        items = os.listdir(sourceDir)
        for item in items:
            if (item.endswith(".dat") or item.endswith(".npy")):
                if not os.access(dataDir+item,os.F_OK):
                    shutil.copy(sourceDir+item,dataDir)
                    self.log.debug("Copied %s to %s" % (item,dataDir))
                else:
                    self.log.debug("File %s already exists in destination %s" % (item,dataDir))
            else:
                self.log.debug("Skipped copying %s, not requried." % (item))
    
    @depends("setup-dirs")    
    def make_massey(self):
        """Make the Massey.FITS file"""
        # Massey.py -> MasseySky.FITS
        from Data.massey import skyab
        # Convert some sky spec from AB to fLambda
        # This sky spectrum is imported from the massey module, and uses the mag conversion functions above.
        skyflam = np.array([skyab[:,0], abmag_to_flambda(skyab[:,1], skyab[:,0])])
        MasseySky = SpectraObject(filename=self.config["Dirs"]["Data"] + "/MasseySky.fits")
        MasseySky.save(skyab.T,"SkyAB")
        MasseySky.save(skyflam,"SkyFL")
        MasseySky.write(clobber=True)
    
    @depends("setup-dirs")        
    def make_hansuchik(self):
        """Make Hansuchik.FITS file"""
        # Hansuchik.py -> HansuchikUVES.FITS
        from Data.Hansuchik import uves_sky
        HansuchikUVES = SpectraObject(filename=self.config["Dirs"]["Data"] + "/UVESSky.fits")
        HansuchikUVES.save(uves_sky.T,"UVESSky")
        HansuchikUVES.write(clobber=True)
    
    @depends("setup-dirs")        
    def make_turnrose(self):
        """Make Turnrose.FITS file"""
        # Turnrose.py -> Turnrose.FITS
        from Data.Turnrose import skyspec
        TurnroseSKY = SpectraObject(filename=self.config["Dirs"]["Data"] + "/TurnroseSKY.fits")
        TurnroseSKY.save(skyspec.T,"TurnroseSKY")
        TurnroseSKY.write(clobber=True)
        
    @depends("setup-dirs")    
    def make_quimby(self):
        """Make Quimby.FITS file"""
        # Quimby.py -> QuimbySky.FITS
        from Data.Quimby import quimby_sky
        QuimbySKY = SpectraObject(filename=self.config["Dirs"]["Data"] + "/QuimbySky.fits")
        QuimbySKY.save(quimby_sky.T,"QuimbySky")
        QuimbySKY.write(clobber=True)
    
    @depends("setup-dirs")    
    def make_atmosphere(self):
        """Make Atmosphere.FITS file"""
        # atmosphere.py -> atmosphere.FITS
        from Data.atmosphere import palextinct
        atmosphereEXT = SpectraObject(filename=self.config["Dirs"]["Data"] + "/atmosphere.fits")
        atmosphereEXT.save(palextinct.T,"Atmosph")
        atmosphereEXT.write(clobber=True)

    @depends("setup-dirs","get-data")    
    def make_palsky(self):
        """Make palsky.FITS file"""
        # palsky_100318.dat -> Palsky.FITS
        palskydata = np.genfromtxt(self.config["Dirs"]["Data"] + "/palsky_100318.dat").T
        palSKY = SpectraObject(filename=self.config["Dirs"]["Data"] + "/PalSky.fits")
        palSKY.save(palskydata,"PalSky")
        palSKY.write(clobber=True)
    
    @depends("setup-dirs")    
    @help("Generate a basic configuration file for this data set.")
    def config_file(self):
        """Basic configrution file."""
        sourceDir = resource_filename(__name__,"/")
        sourceFile = self.config["Configurations"]["This"]
        if not os.access(sourceFile, os.F_OK):
            shutil.copy(sourceDir+sourceFile,".")
            self.log.debug("Copied %s to %s" % (sourceFile,"."))
        else:
            self.log.warning("Not over-writing existing configuration file!")
        
    
def run():
    SA = SetupAgent()
    SA.run()
    
if __name__ == '__main__':    
    run()
        
        
