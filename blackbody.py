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

try:
    
    import AstroObjects
    from AstroObjects.AstroSpectra import SpectraObject
    from AstroObjects.AstroImage import ImageObject
    from AstroObjects.Utilities import BlackBody, get_padding
    
    print("AstroObjects Version %s" % AstroObjects.__version__)
    
except ImportError:
    
    print("Critical Error! Cannot Import AstroObject Model")
    raise
    sys.exit(1)

spectra = SpectraObject()

# Generate a generic BlackBody Spectrum to play with:

# Intial Conditions
temperatures = [4000,5000,6000]
WLstart = 0.37e-6 # In m
WLfinish = 1e-6 # In m
Points = 1e6
WL = np.linspace(WLstart,WLfinish,Points)
PixelSize = 13.5e-3 # In mm
PixelOffset = 2 # In mm, accounts for the distance measurement being from the 0 dispersion point in the optical system
ImageSize = (20,500) # In pixels

for T,i in zip(temperatures,range(len(temperatures))):
    
    # Generate the BlackBody spectrum
    LOG.info("Generating BlackBody Curve at %dK from %.1em to %.1em with %.1e points" % (T,WLstart,WLfinish,Points))
    Flux = BlackBody(WL,T)
    
    # Save the spectrum to our object
    spectra.save(np.array([WL,Flux]),"BlackBody Absolute Units %d" % T)
    
    # Plot the spectra in the normal style
    plt.figure(1)
    plt.plot(WL,Flux,label="BlackBody at %dK" % T)
    plt.xlim([0.35e-6,1.02e-6])
    plt.gca().ticklabel_format(style="sci",scilimits=(3,3))
    plt.xlabel("Wavelength")
    plt.ylabel("J/s")
    plt.legend()
    
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
    # Convert the spectrum to MM using the SPLINE function
    WL,Flux = spectra.data()
    spectra.save(np.array([spline(WL),Flux]),"mm on CCD, spline conversion %d" % T)
    
    # Convert the spectrum to Pixels
    LOG.info("Converting the Spectrum from CCD mm to CCD pixels %d" % T)
    MM,Flux = spectra.data()
    
    # Establish an array of pixels
    Pixels = np.arange(ImageSize[1])
    
    # Find the pixel for each position in MM
    MMPixels = np.digitize((MM+PixelOffset)/PixelSize,Pixels)
    
    LOG.info("Integrating across pixels")
    # Integrate (or simply sum, no functional integration here...)
    PixelFluxes,FPixels = np.histogram(MMPixels,Pixels,weights=Flux)
    
    # Clip out non-zero elements to save
    nonZeroIDX = np.nonzero(PixelFluxes)
    
    # Save the spectra everywhere
    spectra.save(np.array([Pixels[:-1],PixelFluxes]),"Pixels on CCD (all) %d" % T)
    # Save only the data which isn't zero.
    spectra.save(np.array([Pixels[nonZeroIDX],PixelFluxes[nonZeroIDX]]),"Pixels on CCD (nonzero) %d" % T)
    
    # Create a new image
    image = ImageObject()
    
    # The size of the spectra to insert into the image
    SpectraSize = (4,spectra.data()[1].size)
    
    # Generate a flat Image
    data = np.zeros(ImageSize)
    # Insert the spectrum
    data[9:13,spectra.data()[0].astype(int)] = np.resize(spectra.data()[1],SpectraSize)
    
    # Save that generated image
    image.save(data,"GeneartedImage %d" % T)
    
    # Show the image in a field.
    plt.figure(2)
    plt.subplot(len(temperatures),1,i)
    image.show()
    plt.title("Generated Spectrum at %sK" % T)

plt.savefig("BlackbodyImage.png")
plt.figure(1)
plt.savefig("BlackbodySpectrum.png")
plt.show()
LOG.info("Done with Program")

