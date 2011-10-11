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

total = [0.0,0.0]
last = 0.0
start = time.clock()

def timing(string,index):
    """docstring for timing"""
    global total,last,start,stop
    stop = time.clock()
    last = stop - start
    total[index] += last
    LOG.debug(string+" took %fs" % last)
    start = time.clock()

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pyfits
from scipy.interpolate import splprep, splev, splrep, interp1d

timing("Module Import",0)

from AstroObjects.AstroSpectra import FITSSpectra
from AstroObjects.AstroImage import FITSImage
from AstroObjects.Utilities import BlackBody

timing("Custom Module Import",0)

spectra = FITSSpectra()

# Generate a generic BlackBody Spectrum to play with:

T = 5000
WLstart = 0.37e-6
WLfinish = 1e-6
Points = 1e6
LOG.info("Generating BlackBody Curve at %dK from %.1em to %.1em with %.1e points" % (T,WLstart,WLfinish,Points))
WL = np.linspace(WLstart,WLfinish,Points)
Flux = BlackBody(WL,T)

# Save the spectrum to our object
spectra.save(np.array([WL,Flux]),"BlackBody Absolute Units")

timing("Generating BlackBody",0)

# Display the spectrum
LOG.info("Showing the original BlackBody Spectrum")
plt.figure(1)
spectra.showSpectrum()
plt.xlabel("Wavelength")
plt.ylabel("Flux")
plt.title("BlackBody Spectrum at %dK" % T)

timing("Showing BlackBody",1)

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

timing("Generating Conversion Functions",0)

# Demonstrate these functions and thier residuals
LOG.info("Showing Conversion Functions")
shortWL = np.linspace(WLstart,WLfinish,1e4)
plt.figure(2)
plt.plot(shortWL,floor(shortWL),'r.')
plt.plot(shortWL,linear(shortWL),'g.')
plt.plot(shortWL,spline(shortWL),'b.')
plt.title("Wavelength Conversion Functions")
plt.xlabel("Wavelength (m)")
plt.ylabel("Distance Along CCD (microns)")
plt.gca().ticklabel_format(style="sci",scilimits=(3,3))

timing("Showing Conversion Functions",1)

LOG.info("Showing Wavelength Residuals")
LOG.info("Note: This doesn't make sense for floor and linear interpolation, as they will always be exact at given data points!")
givenWL = WLtoMM[:,0]
givenMM = WLtoMM[:,1]
plt.figure(3)
plt.plot(givenWL,floor(givenWL*1e-6)-givenMM,'r.')
plt.plot(givenWL,linear(givenWL*1e-6)-givenMM,'g.')
plt.plot(givenWL,spline(givenWL*1e-6)-givenMM,'b.')
plt.title("Wavelength Conversion Functions Residuals")
plt.xlabel("Wavelength (m)")
plt.ylabel("Residual Distance Along CCD (mm)")
plt.gca().ticklabel_format(style="sci",scilimits=(3,3))

timing("Showing Conversion Function Residuals",1)

LOG.info("Converting the Spectrum from Wavelength to CCD mm")
# Convert the spectrum to MM
WL,Flux = spectra.data()
spectra.save(np.array([spline(WL),Flux]),"mm on CCD, spline conversion")

timing("Converting to MM via Spline",0)

LOG.info("Showing the converted Spectrum")
plt.figure(4)
spectra.showSpectrum()
plt.xlabel("CCD Distance (mm)")
plt.ylabel("Flux")
plt.title("BlackBody Spectrum at %dK" % T)

timing("Showing Converted Spectrum in mm",1)

# Convert the spectrum to Pixels
LOG.info("Converting the Spectrum from CCD mm to CCD pixels")
MM,Flux = spectra.data()

PixelSize = 13.5e-3
PixelOffset = 2
ImageSize = (20,500)

# Establish an array of pixels
Pixels = np.arange(ImageSize[1])

# Find the pixel for each position in MM
MMPixels = np.digitize((MM+PixelOffset)/PixelSize,Pixels)

timing("Genearting Pixels for Conversion",0)


LOG.info("Integrating across pixels")
# Integrate (or simply sum, no functional integration here...)
PixelFluxes,FPixels = np.histogram(MMPixels,Pixels,weights=Flux)

timing("Integrating over pixels for conversion",0)

# Clip out non-zero elements to save
nonZeroIDX = np.nonzero(PixelFluxes)

spectra.save(np.array([Pixels,PixelFluxes]),"Pixels on CCD (all)")
spectra.save(np.array([Pixels[nonZeroIDX],PixelFluxes[nonZeroIDX]]),"Pixels on CCD (nonzero)")

timing("Saving spectrum in pixels",0)


LOG.info("Showing the pixelated spectrum")
plt.figure(5)
spectra.showSpectrum()
plt.xlabel("CCD Pixels")
plt.ylabel("Flux")
plt.title("BlackBody Spectrum at %dK" % T)
plt.gca().ticklabel_format(style="plain",axis="x")

timing("Showing Spectrum in pixels",1)

image = FITSImage()
SpectraSize = (4,spectra.data()[1].size)

data = np.zeros(ImageSize)
data[9:13,spectra.data()[0].astype(int)] = np.resize(spectra.data()[1],SpectraSize)

image.save(data,"GeneartedImage")
timing("Generating Image",0)

plt.figure(6)
image.show()
plt.title("Generated Spectrum in Place")
timing("Showing Image",1)

LOG.info("Data functions took %f to execute" % total[0])
LOG.info("Plotting functions took %f to execute" % total[1])
plt.show()


