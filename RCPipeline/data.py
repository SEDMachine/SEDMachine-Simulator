#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  data.py
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

from AstroObject.AstroImage import ImageObject
from AstroObject.AstroObjectLogging import logging
exptime = 1200

images = ImageObject()
images["BiasA"] = np.random.poisson(10,(1000,1000)).astype(np.uint16)
images["BiasB"] = np.random.poisson(10,(1000,1000)).astype(np.uint16)
images["BiasC"] = np.random.poisson(10,(1000,1000)).astype(np.uint16)
images["DarkA"] = np.random.poisson(10,(1000,1000)).astype(np.uint16)
images["DarkB"] = np.random.poisson(10,(1000,1000)).astype(np.uint16)
images["DarkC"] = np.random.poisson(10,(1000,1000)).astype(np.uint16)
for image in ["DarkA","DarkB","DarkC"]:
    images[image].header.update('exptime',exptime)
images["FlatA"] = np.ones((1000,1000)).astype(np.uint16)
images["FlatB"] = np.ones((1000,1000)).astype(np.uint16)
images["FlatC"] = np.ones((1000,1000)).astype(np.uint16)
images["DataA"] = np.ones((1000,1000)).astype(np.uint16)
images["DataB"] = np.ones((1000,1000)).astype(np.uint16)
images["DataC"] = np.ones((1000,1000)).astype(np.uint16)
for image in ["DataA","DataB","DataC"]:
    images[image].header.update('exptime',exptime)


images.write("Bias.fits",states=["BiasA","BiasB","BiasC"],clobber=True)
images.write("Dark.fits",states=["DarkA","DarkB","DarkC"],clobber=True)
images.write("Flat.fits",states=["FlatA","FlatB","FlatC"],clobber=True)
images.write("Data.fits",states=["DataA","DataB","DataC"],clobber=True)