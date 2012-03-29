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

from AstroObject.AstroSimulator import Simulator
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
        self.config.update({"Default":"*all"})
        self._register_stages()
        
    def _register_stages(self):
        """Register all of the required stages"""
        self.registerStage(self.EnsureDirectories,"setup-dirs",description="Setting up directories",include=True,help="Make required Directories")
        self.registerStage(self.GetStaticData,"get-data",description="Collecting static (.dat) files",include=True,dependencies=["setup-dirs"],help="Populate data directory with included static data")
        self.registerStage(self.MakeMassey,"make-massey",dependencies=["setup-dirs"])
        self.registerStage(self.MakeHansuchik,"make-hansuchik",dependencies=["setup-dirs"])
        self.registerStage(self.MakeTurnrose,"make-turnrose",dependencies=["setup-dirs"])
        self.registerStage(self.MakeQuimby,"make-quimby",dependencies=["setup-dirs"])
        self.registerStage(self.MakeAtmosphere,"make-atmosphere",dependencies=["setup-dirs"])
        self.registerStage(self.MakePalSky,"make-palsky",dependencies=["setup-dirs","get-data"])
        self.registerStage(None,"make-specs",dependencies=["make-massey","make-hansuchik","make-turnrose","make-quimby","make-atmosphere","make-palsky"],include=True,help="Setup spectral FITS files")
        
    def EnsureDirectories(self):
        """Ensure that the required directories exist."""
        for key,directory in self.config["Dirs"].iteritems():
            try:
                os.mkdir(directory)
            except os.error, e:
                self.log.debug("Directory %s already exists." % directory)
        
    def GetStaticData(self):
        """Copy static data files to the data directory."""
        dataDir = self.config["Dirs"]["Data"] + "/"
        sourceDir = resource_filename(__name__,"Data/")
        items = os.listdir(sourceDir)
        for item in items:
            if item.endswith(".dat") and not os.access(dataDir+item,os.F_OK):
                shutil.copy(sourceDir+item,dataDir)
                self.log.debug("Copied %s to %s" % (item,dataDir))
            else:
                self.log.debug("Skipped copying %s" % (item))
        
    def MakeMassey(self):
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
        
    def MakeHansuchik(self):
        """Make Hansuchik.FITS file"""
        # Hansuchik.py -> HansuchikUVES.FITS
        from Data.Hansuchik import uves_sky
        HansuchikUVES = SpectraObject(filename=self.config["Dirs"]["Data"] + "/UVESSky.fits")
        HansuchikUVES.save(uves_sky.T,"UVESSky")
        HansuchikUVES.write(clobber=True)
        
    def MakeTurnrose(self):
        """Make Turnrose.FITS file"""
        # Turnrose.py -> Turnrose.FITS
        from Data.Turnrose import skyspec
        TurnroseSKY = SpectraObject(filename=self.config["Dirs"]["Data"] + "/TurnroseSKY.fits")
        TurnroseSKY.save(skyspec.T,"TurnroseSKY")
        TurnroseSKY.write(clobber=True)
        
    def MakeQuimby(self):
        """Make Quimby.FITS file"""
        # Quimby.py -> QuimbySky.FITS
        from Data.Quimby import quimby_sky
        QuimbySKY = SpectraObject(filename=self.config["Dirs"]["Data"] + "/QuimbySky.fits")
        QuimbySKY.save(quimby_sky.T,"QuimbySky")
        QuimbySKY.write(clobber=True)
    
    def MakeAtmosphere(self):
        """Make Atmosphere.FITS file"""
        # atmosphere.py -> atmosphere.FITS
        from Data.atmosphere import palextinct
        atmosphereEXT = SpectraObject(filename=self.config["Dirs"]["Data"] + "/atmosphere.fits")
        atmosphereEXT.save(palextinct.T,"Atmosph")
        atmosphereEXT.write(clobber=True)

    def MakePalSky(self):
        """Make palsky.FITS file"""
        # palsky_100318.dat -> Palsky.FITS
        palskydata = np.genfromtxt(self.config["Dirs"]["Data"] + "/palsky_100318.dat").T
        palSKY = SpectraObject(filename=self.config["Dirs"]["Data"] + "/PalSky.fits")
        palSKY.save(palskydata,"PalSky")
        palSKY.write(clobber=True)
    
def run():
    SA = SetupAgent()
    SA.run()
    
if __name__ == '__main__':    
    run()
        
        
