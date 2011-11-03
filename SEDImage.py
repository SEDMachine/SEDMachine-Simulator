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
import scipy.signal

from scipy.interpolate import splprep, splev, splrep, interp1d

LOG = logging.getLogger("Basic SED Image")

import AstroObject
from AstroObject.AstroSpectra import SpectraObject
from AstroObject.AstroImage import ImageObject,ImageFrame
from AstroObject.AnalyticSpectra import BlackBodySpectrum, AnalyticSpectrum, FlatSpectrum
from AstroObject.Utilities import *

from Utilities import *

class SEDSystem(object):
    """A Model for the SED System which can be used by SEDImage to place spectra etc."""
    def __init__(self):
        super(SEDSystem, self).__init__()
        # Some data about this system
        self.spectrumFWHM = 2/1.26 # In pixels
        self.widthmm = 35.0
        self.pixunit = 1 / (0.0135/1) #Pixels per mm
        self.npix = np.round(self.widthmm * self.pixunit, 0) #Number of Pixels
        self.deltmm = self.widthmm / self.npix #Spacing, in mm, between pixels
        self.PADDING = 5
        
    def gauss_kern(self,size, sizey=None):
        """ Returns a normalized 2D gauss kernel array for convolutions """
        size = int(size)
        if not sizey:
            sizey = size
        else:
            sizey = np.int(sizey)
        x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
        g = np.exp(-(x**2/np.float(size)+y**2/np.float(sizey)))
        return g / g.sum()
    
    def circular_kern(self,radius):
        """Returns a circular kernel of the given radius"""
        radius = int(radius)
        
        x,y = np.mgrid[-radius:radius+1, -radius:radius+1]
        d = np.sqrt(x**2.0 + y**2.0)
        tv = d <= radius
        return tv.astype(np.float)
        
    def blur_image(self,im, n, ny=None) :
        """ blurs the image by convolving with a gaussian kernel of typical
            size n. The optional keyword argument ny allows for a different
            size in the y direction.
        """
        g = self.gauss_kern(n, sizey=ny)
        improc = sp.signal.convolve(im,g, mode='same')
        return(improc)
        
    def LoadDispersionData(self,dispspec):
        """This loads the dispersion data. Dispersion data is provided in mm along the spectrum and wavelength. This method saves both the raw data (in :attr:`self.MM` and :attr:`self.WL`) as well as functions which can be used to convert across the dispersion. Functions by default take a mm position and return the wavelength at this position."""
        WLtoMM = np.genfromtxt(dispspec).T
        WL = WLtoMM[0]
        MM = WLtoMM[1]        
        MM -= np.min(MM) #Zero offset the MM data
        self.MM,self.WL = MM,WL
        
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
        self.spvars = splrep(MM,WL)
        self.spline = lambda x: splev(x,self.spvars)
        
    def LoadLensletData(self,laspec):
        """docstring for Load Lenselt Data"""
        # Load Lenslet Specification File
        ix, p1, p2, lams, xs, ys = np.genfromtxt(laspec,skip_header=1).T
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
        xs += (self.widthmm/2)
        ys += (self.widthmm/2)
        
        # Find the xs and ys that are not within 0.1 mm of the edge of the detector...
        ok = (xs > 0.1) & (xs < self.widthmm-0.1) & (ys > 0.1) & (ys < self.widthmm-0.1)
        ix, p1, p2, lams, xs, ys = ix[ok], p1[ok], p2[ok], lams[ok], xs[ok], ys[ok]
        # We remove these positions because we don't want to generate spectra for them.
        
        # This simply generates a list of all of the lenslets
        self.lenslets = np.unique(ix)
        
        # Convert the xs and ys to pixel positions
        xpix = np.round(xs * self.pixunit,0).astype(np.int)
        ypix = np.round(ys * self.pixunit,0).astype(np.int)
        
        # Determine the center of the whole system by finding the x position that is closest to 0,0 in pupil position
        cntix = np.argmin(p1**2 + p2**2)
        self.center = (xs[cntix] * self.pixunit, ys[cntix] * self.pixunit)
        
        self.ix, self.p1, self.p2, self.lams, self.xs, self.ys = ix, p1, p2, lams, xs, ys
        self.xpix, self.ypix = xpix, ypix
        
    def loadOpticsData(self,laspec,dispspec):
        """Loads an optical conversion based on the lenslet array spec and dispersion spec files provided"""
        
        self.LoadDispersionData(dispspec)
        self.LoadLensletData(laspec)
        
        
        
    def get_dense_image(self,lenslet,spectrum,model,density):
        """This function returns a dense image array, with flux placed into single pixels."""
        points, intpoints, wl, deltawl, density = model.get_wavelengths(lenslet,density)
        radiance = spectrum(wl*1e-6) * 1e-6
        flux = radiance[1,:-1] * deltawl

        x,y = intpoints[:-1].T.astype(np.int)

        xint,yint = points.T.astype(np.int)

        x -= np.min(x)
        y -= np.min(y)

        xdist = np.max(x)-np.min(x)
        ydist = np.max(y)-np.min(y)

        xdist += (5 - xdist % 5)
        ydist += (5 - ydist % 5)

        x += self.PADDING * density
        y += self.PADDING * density

        img = np.zeros((xdist+2*self.PADDING*density,ydist+2*self.PADDING*density))

        img[x,y] = flux

        corner = np.array([ xint[np.argmax(x)], yint[np.argmin(y)]]) - np.array([self.PADDING,self.PADDING])

        return img, corner
    
    def get_wavelengths(self,lenslet_num,density=1):
        """Returns an array of wavelengths for a given lenslet number"""
        # This method will be the critical slow point for the whole system...
        # It should be re-written to do the conversions in a better way, but
        # I'm not sure how to do that right now.
        
        
        # First, we take only data points which apply to this lenslet
        use = lenslet_num == self.ix
        
        # Each lenslet should have three points. It might not though, because we clipped some lenselet points that were too close to the edge. 
        # We really should move this logic up higher.
        if len(self.xpix[use]) < 3: 
            raise SEDLimits
        
        if np.any(self.xpix[use] == 0):
            raise SEDLimits
        
        if np.any(self.ypix[use] == 0):
            raise SEDLimits
        
        # We ignore anything that has an x-dispersion across more than 30 pixels.
        if np.any(np.abs(np.diff(self.xpix[use])) > 30):
            raise SEDLimits

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
        
        if distance == 0:
            raise SEDLimits
        
        # Find the length in units of (int) pixels
        npix = (distance * self.pixunit).astype(np.int) * density

        # Create a data array one hundred times as dense as the number of pixels
        #   This is the super dense array which will use the above interpolation
        dummy_lam = np.linspace(np.min(self.lams[use]),np.max(self.lams[use]),npix*100)

        # Interpolate along our really dense set of wavelengths to find all possible
        # illuminated pixel positions in this spectrum
        dummy_pts = np.array([fx(dummy_lam),fy(dummy_lam)])
        
        if density == 1:
            dummy_pts = dummy_pts.astype(np.int)
        else:
            dummy_pts = np.round(dummy_pts * density) / density
        
        # Remove any duplicate points. This does not do so in order, so we must
        # sort the array of unique points afterwards...
        unique_pts = np.array([np.array(x) for x in set(tuple(x) for x in dummy_pts.T)])
        
        # An array of distances to the origin of this spectrum, can be used to find wavelength
        # of light at each distance
        distance = np.array([np.sqrt(np.sum((x/self.pixunit)-start)**2) for x in unique_pts])

        sorted_idx = np.argsort(distance)

        distance = distance[sorted_idx]
        points = unique_pts[sorted_idx]
        intpoints = points * density
        
        wl = self.spline(distance)
        if density==1:
            return points,wl,np.diff(wl)
        else:
            return points,intpoints,wl,np.diff(wl),density
    

class SEDImage(ImageObject):
    """Representation of an SEDMachine Image, with methods for placing spectra"""
    
    def get_sub_image(self,lenslet,spectrum,model,density,stdev):
        """docstring for get_sub_image"""
        img,corner = model.get_dense_image(lenslet,spectrum,model,density)
        img2 = model.blur_image(img,stdev*density)
        small = bin(img2,density).astype(np.int16)
        
        return small,corner
    
    def place_sed(self,lenslet,spectrum,label,model,density,stdev):
        """docstring for place_sed"""
        small, corner = self.get_sub_image(lenslet,spectrum,model,density,stdev)
        
        self.place(small,corner,label)
        
    def place(self,img,corner,label):
        """Place the given AstroObject.AnalyticSpectra.AnalyticSpectrum onto the SEDMachine Image"""
        
        xstart = corner[0]
        xend = xstart + img.shape[0]
        ystart = corner[1]
        yend = ystart + img.shape[1]
        data = self.data()
        
        if data.shape[0] < xend or data.shape[1] < yend:
            raise SEDLimits
        
        if xstart < 0 or ystart < 0:
            raise SEDLimits
        
        if xend < 0 or yend < 0:
            raise SEDLimits
        
        data[xstart:xend,ystart:yend] += img
        
        
    def generate_blank(self,shape):
        """Generates a blank SEDMachine Image"""
        self.save(np.zeros(shape).astype(np.int16),"Blank")
        
        
    def crop(self,x,y,xsize,ysize=None):
        """Crops the provided image to twice the specified size, centered around the x and y coordinates provided."""
        if not ysize:
            ysize = xsize
        cropped = self.states[self.statename].data[x-xsize:x+xsize,y-ysize:y+ysize]
        LOG.debug("Cropped and Saved Image")
        self.save(cropped,"Cropped")
        

class SEDLimits(Exception):
    """docstring for SEDLimits"""
    pass
        
        