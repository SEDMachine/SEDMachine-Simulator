#!/usr/bin/env python
# 
#  plot.py
#  SEDMachine
#  
#  Created by Alexander Rudy on 2011-10-12.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 

import time,logging

LOG = logging.getLogger("Basic Operations")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pyfits
from scipy.interpolate import splprep, splev, splrep, interp1d

from AstroObjects.AstroSpectra import FITSSpectra
from AstroObjects.AstroImage import FITSImage
from AstroObjects.Utilities import BlackBody

spectra = FITSSpectra()

# Generate a generic BlackBody Spectrum to play with:

temperatures = [4000,5000,6000]

for T,i in zip(temperatures,range(len(temperatures))):
    WLstart = 0.37e-6
    WLfinish = 1e-6
    Points = 1e6
    LOG.info("Generating BlackBody Curve at %dK from %.1em to %.1em with %.1e points" % (T,WLstart,WLfinish,Points))
    WL = np.linspace(WLstart,WLfinish,Points)
    Flux = BlackBody(WL,T)
    
    # Save the spectrum to our object
    spectra.save(np.array([WL,Flux]),"BlackBody Absolute Units %d" % T)
    
    plt.figure(1)
    plt.plot(WL,Flux)
    plt.gca().ticklabel_format(style="sci",scilimits=(3,3))
    plt.xlabel("Wavelength")
    plt.ylabel("Joules")
    
    LOG.info("Constructing Conversion from WL to mm along CCD...")
    # Read the conversion file
    WLtoMM = np.genfromtxt("Rtest.dat")
    
    # Define a basic Floor function interpolation
    floor = lambda x: np.vectorize(lambda ax: WLtoMM[(np.abs(ax-(WLtoMM[:,0]*1e-6))).argmin(),1])(x)
    # Define a piecewise linear interpolation
    linear = interp1d(WLtoMM[:,0]*1e-6,WLtoMM[:,1])
    # Generate the SPLINE function:
    spvars = splrep(WLtoMM[:,0]*1e-6,WLtoMM[:,1])
    spline = lambda x: splev(x,spvars)
    
    LOG.info("Converting the Spectrum from Wavelength to CCD mm")
    # Convert the spectrum to MM
    WL,Flux = spectra.data()
    spectra.save(np.array([spline(WL),Flux]),"mm on CCD, spline conversion %d" % T)
    
    # Convert the spectrum to Pixels
    LOG.info("Converting the Spectrum from CCD mm to CCD pixels %d" % T)
    MM,Flux = spectra.data()
    
    PixelSize = 13.5e-3
    PixelOffset = 2
    ImageSize = (20,500)
    
    # Establish an array of pixels
    Pixels = np.arange(ImageSize[1])
    
    # Find the pixel for each position in MM
    MMPixels = np.digitize((MM+PixelOffset)/PixelSize,Pixels)
    
    LOG.info("Integrating across pixels")
    # Integrate (or simply sum, no functional integration here...)
    PixelFluxes,FPixels = np.histogram(MMPixels,Pixels,weights=Flux)
    
    # Clip out non-zero elements to save
    nonZeroIDX = np.nonzero(PixelFluxes)
    
    spectra.save(np.array([Pixels[:-1],PixelFluxes]),"Pixels on CCD (all) %d" % T)
    spectra.save(np.array([Pixels[nonZeroIDX],PixelFluxes[nonZeroIDX]]),"Pixels on CCD (nonzero) %d" % T)
    
    image = FITSImage()
    SpectraSize = (4,spectra.data()[1].size)
    
    data = np.zeros(ImageSize)
    data[9:13,spectra.data()[0].astype(int)] = np.resize(spectra.data()[1],SpectraSize)
    
    image.save(data,"GeneartedImage %d" % T)
    
    plt.figure(2)
    plt.subplot(len(temperatures),1,i)
    image.show()
    plt.title("Generated Spectrum at %sK" % T)
plt.show()
LOG.info("Done with Program")

