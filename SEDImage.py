#!/usr/bin/env python
# 
#  image.py
#  SEDMachine
#  
#  Created by Alexander Rudy on 2011-10-19.
#  Copyright 2011 Alexander Rudy. All rights reserved.
#  This file generates a simple FITS file with blackbody spectra at the position of each spectrum specified for SED Machine.


import logging
import arpytools.progressbar

import numpy as np
import pyfits as pf
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.interpolate import splprep, splev, splrep, interp1d

LOG = logging.getLogger("Basic SED Image")

import AstroObject
from AstroObject.AstroSpectra import SpectraObject
from AstroObject.AstroImage import ImageObject,ImageFrame
from AstroObject.AnalyticSpectra import BlackBodySpectrum, AnalyticSpectrum
from AstroObject.Utilities import BlackBody, get_padding

from AOFuture import expandLim, disable_Console, enable_Console



class SEDImage(ImageObject):
    """Representation of an SEDMachine Image, with methods for placing spectra"""
    
    
    def loadOpticsData(self,laspec,dispspec):
        """Loads an optical conversion based on the lenslet array spec and dispersion spec files provided"""
        # Some data
        
        self.widthmm = 35.0 #Width of the Focal Plane
        self.pixunit = 1 / (0.0135/1) #Pixels per mm
        self.npix = np.round(self.widthmm * self.pixunit, 0) #Number of Pixels
        self.deltmm = self.widthmm / self.npix #Spacing, in mm, between pixels
        
        
        # Load Lenslet Specification File
        ix, p1, p2, lams, xs, ys = np.genfromtxt(laspec,skip_header=1).T

        # Correctly Type Lenslet Specification Data
        
        ix = ix.astype(np.int) #Indicies should always be integers
        
        # Center xs and ys on detector with a corner at 0,0
        xs += self.widthmm/2
        ys += self.widthmm/2
        
        # Find the xs and ys that are not within 0.1 mm of the edge of the detector...
        ok = (xs > 0.1) & (xs < self.widthmm-0.1) & (ys > 0.1) & (ys < self.widthmm-0.1)
        ix, p1, p2, lams, xs, ys = ix[ok], p1[ok], p2[ok], lams[ok], xs[ok], ys[ok]
        # Pupils are 2.4 px in diameter
        
        self.lenslets = np.unique(ix)
        
        xpix = np.round(xs * self.pixunit,0).astype(np.int)
        ypix = np.round(ys * self.pixunit,0).astype(np.int)
        
        cntix = np.argmin(p1**2 + p2**2)
        self.center = (xs[cntix] * self.pixunit, ys[cntix] * self.pixunit)
        
        WLtoMM = np.genfromtxt("Rtest.dat").T
        
        WL = WLtoMM[0]
        MM = WLtoMM[1]
        
        MM -= np.min(MM)
        
        # A note about functions:
        #
        # All of these functions will return the wavelength of a position, in milimeters,
        # from the origin of the spectrum. We define the origin of the spectrum as the
        # postion where the 0.37 micron wavelength falls.
        
        # Define a basic Floor function interpolation
        self.floor = lambda x: np.vectorize(lambda ax: MM[(np.abs(ax-WL)).argmin()])(x)
        # Define a piecewise linear interpolation
        self.linear = interp1d(WL,MM)
        # Generate the SPLINE function:
        spvars = splrep(MM,WL)
        self.spline = lambda x: splev(x,spvars)
        
        
        self.MM,self.WL = MM,WL
        self.ix, self.p1, self.p2, self.lams, self.xs, self.ys = ix, p1, p2, lams, xs, ys
        self.xpix, self.ypix = xpix,ypix
    
    def get_wavelengths(self,lenslet_num):
        """Returns an array of wavelengths for a given lenslet number"""
        # This method will be the critical slow point for the whole system...
        # It should be re-written to do the conversions in a better way, but
        # I'm not sure how to do that right now.
        
        
        # First, we take only data points which apply to this lenslet
        use = lenslet_num == self.ix
        
        if len(self.xpix[use]) < 3: 
            raise Exception

        if np.any(self.xpix[use] == 0):
            raise Exception

        if np.any(self.ypix[use] == 0):
            raise Exception

        if np.any(np.abs(np.diff(self.xpix[use])) > 30):
            raise Exception
        
        # Interpolation to convert from wavelength to pixels.
        #   The accuracy of this interpolation is not important.
        #   Rather, it is used to find the pixels where the light will fall
        #   and is fed an array that is very dense, used on this dense interpolation
        #   and then binned back onto pixels. Thus it will be used to get a list
        #   of all illuminated pixels.
        fx = np.poly1d(np.polyfit(self.lams[use], self.xpix[use], 2))
        fy = np.poly1d(np.polyfit(self.lams[use], self.ypix[use], 2))
        
        # Find the starting and ending position of the spectra
        startix = np.argmin(self.lams[use])
        endix = np.argmax(self.lams[use])
        start = np.array([self.xs[use][startix],self.ys[use][startix]])
        end = np.array([self.xs[use][endix],self.ys[use][endix]])
        
        # Get the total length of the spectra
        distance = np.sqrt(np.sum(end-start)**2)
        
        # Find the length in units of (int) pixels
        npix = (distance * self.pixunit).astype(np.int)
        
        # Create a data array one hundred times as dense as the number of pixels
        #   This is the super dense array which will use the above interpolation
        dummy_lam = np.linspace(np.min(self.lams[use]),np.max(self.lams[use]),npix*10)
        
        # Interpolate along our really dense set of wavelengths to find all possible
        # illuminated pixel positions in this spectrum
        dummy_pts = np.array([fx(dummy_lam).astype(np.int),fy(dummy_lam).astype(np.int)])
        
        # Remove any duplicate points. This does not do so in order, so we must
        # sort the array of unique points afterwards...
        unique_pts = np.array([np.array(x).T for x in set(tuple(x) for x in dummy_pts.T)]).T
        
        # An array of distances to the origin of this spectrum, can be used to find wavelength
        # of light at each distance
        distance = np.array([np.sqrt(np.sum((x/self.pixunit)-start)**2) for x in unique_pts.T])
        
        sorted_idx = np.argsort(distance)
        
        distance = distance[sorted_idx]
        points = unique_pts[:,sorted_idx]
        
        wl = self.spline(distance)
        
        return points,wl,np.diff(wl)
    
    def place(self,lenslet,spectrum,label):
        """Place the given AstroObject.AnalyticSpectra.AnalyticSpectrum onto the SEDMachine Image"""
        image = self.data()
        
        points, wl, deltawl = self.get_wavelengths(lenslet)
        
        radiance = spectrum(wl*1e-6)
                
        flux = radiance[1,:-1] * deltawl
        
        x,y = points[:,:-1]
        
        gain = 1e-10
        
        flux = (flux  * gain).astype(np.int)
        
        image[x,y] = flux
        
        self.save(image,label)
        
    def generate_blank(self):
        """Generates a blank SEDMachine Image"""
        self.save(np.zeros((self.npix, self.npix)).astype(np.int),"Blank")
        
        
    def crop(self,x,y,xsize,ysize=None):
        """Crops the provided image to twice the specified size, centered around the x and y coordinates provided."""
        if not ysize:
            ysize = xsize
        cropped = self.states[self.statename].data[x-xsize:x+xsize,y-ysize:y+ysize]
        LOG.debug("Cropped and Saved Image")
        self.save(cropped,"Cropped")


        