# 
#  Spectra.py
#  SEDMachine
#  
#  Created by Alexander Rudy on 2011-10-04.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pyfits
import logging
from Utilities import *

LOG = logging.getLogger(__name__)

class Spectrum(object):
    """A simple spectrum object"""
    
    xvals = ["Wavelength","Frequency","Energy","Pixels","Distance"]
    yvals = ["Flux","Counts"]
    
    def __init__(self,spectrum=None,x="Wavelength",y="Flux",label="Raw Spectrum"):
        super(Spectrum, self).__init__()
        self.spectrum = None
        if spectrum != None:
            self.save(spectrum,x,y,label)
    
    def save(self,spectrum,x,y,label):
        """Saves the spectrum to the spectrum object"""
        if self.spectrum != None:
            LOG.warn("Overwriting Spectrum")
        if type(spectrum) != np.ndarray:
            raise ValueError("Spectrum must be numpy.ndarray, not %s" % type(spectrum))
        if spectrum.shape[0] != 2:
            raise ValueError("Spectrum object only accepts 1-D spectra")
        self.__check_axis(x,y)
        self.spectrum = spectrum
        self.x = x
        self.y = y
        self.label = label
        LOG.info("Saved spectrum data \"%s\" with shape %s." % (label,spectrum.shape))
    
    def generate_ramp(self,start,end,spacing,label="Ramp Spectrum",x="Wavelength",y="Flux"):
        """Generates and saves a linear ramp spectrum with given boundaries"""
        counts = np.linspace(start,end,spacing)
        spectra = np.array([counts,counts])
        LOG.info("Generated Ramp Spectrum from %(bottom)s to %(top)s labeled \"%(label)s\"." % {"bottom":start,"top":end,"label":label})
        self.save(spectra,x,y,label)
        
    def __check_axis(self,x=None,y=None):
        """Checks the given items for axis validity"""
        if x not in self.xvals and x != None:
            raise ValueError("X-Axis not understood: %s \nAllowd Values:%s" % (x,self.xvals))
        if y not in self.yvals and y != None:
            raise ValueError("Y-Axis not understood: %s \nAllowd Values:%s" % (y,self.yvals))
        return
    
    def show_spectrum(self):
        """Applies the necessary plotting items to show a spectrum"""
        x = self.spectrum[0]
        value = self.spectrum[1]
        plt.plot(x,value,'k-')
        set_axis_padding((x,value))
        plt.title("Spectrum: %s" % self.label)
        plt.xlabel(self.x)
        plt.ylabel(self.y)
    
    def load_conversion(self,filename,origin,result,axis="x"):
        """Loads a conversion file of equivalent points along an axis to convert between two different axis types."""
        if axis == "x":
            self.__check_axis(x=origin)
            self.__check_axis(x=result)
        elif axis == "y":
            self.__check_axis(y=origin)
            self.__check_axis(y=result)
        else:
            raise ValueError("Axis Invalid: %s, should be \"x\" or \"y\"!")
        
        conversionFile = open(filename,"r")
    

    
            