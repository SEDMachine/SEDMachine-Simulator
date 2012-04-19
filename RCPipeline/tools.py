#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  tools.py
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2012-04-19.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

# Python Imports
import shutil
import os
import collections

# Numpy Imports
import numpy as np

# Package Resources Imports
from pkg_resources import resource_filename

# PyRAF Imports
from pyraf import iraf
from iraf import imred, ccdred

from AstroObject.AstroSimulator import Simulator
from AstroObject.AstroSimulator import (
    optional,
    description,
    include,
    replaces,
    depends,
    excepts,
    collect,
    ignore,
    help
)
from AstroObject.AstroObjectLogging import logging

class SetupAgent(Simulator):
    """An agent for setting up directories for simulations."""
    def __init__(self):
        super(SetupAgent, self).__init__(commandLine=True)
        self.config.load(resource_filename("RCPipeline","Defaults.yaml"))
        self.config.setFile("Main")
        self.config.load()
        self.config.update({"Default":"*all"})
        self.collect()
        self.registerStage(None,"test",dependencies=["make-test-data","config-file"],help="Set up the test pipeline run",description="Test Mode Setup")
         
    @include
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
    
    
    @help("Generate test Data")
    def make_test_data(self):
        """Making test data files"""
        import data
    
    @help("Generate a Test configuration file for the test data")
    def config_file(self):
        """Basic configrution file."""
        sourceDir = resource_filename(__name__,"/")
        sourceFile = self.config["Configurations"]["Test"]
        destFileName = self.config["Configurations"]["Main"]
        if not os.access(destFileName, os.F_OK):
            shutil.copy(sourceDir+sourceFile,destFileName)
            self.log.debug("Copied %s to %s" % (sourceFile,"./"+destFileName))
        else:
            self.log.warning("Not over-writing existing configuration file!")
    
def main():
    Agent = SetupAgent()
    Agent.run()
