#!/usr/bin/env python
# 
#  SED.py
#  Simulation Software
#  
#  Created by Alexander Rudy on 2011-11-04.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 


import numpy as np
import pyfits as pf
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.signal
import scipy.interpolate

import AstroObject
from AstroObject.AstroSpectra import SpectraObject
from AstroObject.AstroImage import ImageObject,ImageFrame
from AstroObject.AnalyticSpectra import BlackBodySpectrum, AnalyticSpectrum, FlatSpectrum
from AstroObject.Utilities import *

from Utilities import *

class SEDLimits(Exception):
    """A Basic Error-Differentiation Class"""
    pass

class Model(ImageObject):
    """This object is model"""
    def __init__(self):
        super(Model, self).__init__()
        
        # Configuration Variables for The System
        
        self.convert = {}
        self.ccd_size = {}
        self.image_sz = {}
        self.tel_radi = {}
        self.tel_obsc = {}
        self.psf_stdv = {}
        
        
        self.convert["pxtomm"] = 0.0135
        self.convert["mmtopx"] = 1.0 / 0.0135
        
        # CCD / Image Plane Information
        
        self.ccd_size["px"] = 2048 #pixels
        self.ccd_size["mm"] = self.ccd_size["px"] * self.convert["pxtomm"]
        self.image_sz["mm"] = 35.0 #mm
        self.image_sz["px"] = np.round( self.image_sz["mm"] * self.convert["mmtopx"] , 0 )
        
        # Telescope Information
        
        self.tel_radi["px"] = 2.4 / 2.0
        self.tel_radi["mm"] = self.tel_radi["px"] * self.convert["pxtomm"]
        self.tel_obsc["px"] = 0.4 / 2.0
        self.tel_obsc["mm"] = self.tel_obsc["px"] * self.convert["pxtomm"]
        
        # PSF Information
        # For a gaussian PSF
        self.psf_stdv["px"] = 1.0
        
        # Image Generation Density
        self.density = 5
        self.padding = 5
        
        self.gain = 1e-6
        
        self.files = {}
        self.files["lenslets"] = "Data/lenslet_xy_13-10-2011.txt"
        self.files["dispersion"] = "Data/dispersion_12-10-2011.txt"
        
    def setup(self):
        """Function handles the setup of simulation information"""
        
        # Telescope Image Setup
        self.TELIMG = self.get_tel_kern()
        
        # PSF Setup
        self.PSFIMG = self.get_psf_kern()
        
        # Preconvolved System
        self.FINIMG = sp.signal.convolve(self.PSFIMG,self.TELIMG,mode='same')
        
        # Get layout data from files
        self.loadOpticsData(self.files["lenslets"],self.files["dispersion"])
        
        self.generate_blank()
        
    def get_tel_kern(self):
        """Returns the telescope kernel"""
        TELIMG = self.circle_kern( self.tel_radi["px"] * self.density )
        center = self.circle_kern( self.tel_obsc["px"] * self.density , self.tel_radi["px"] * self.density , False )
        TELIMG -= center
        TELIMG = TELIMG / np.sum(TELIMG)
        
        return TELIMG
        
    def get_psf_kern(self):
        """Returns the PSF Kernel"""
        return self.gauss_kern( (self.psf_stdv["px"] * self.density) )
        
    def get_blank_img(self):
        """Returns an image of the correct size to use for spectrum placement"""
        return np.zeros((self.image_sz["px"],self.image_sz["px"]))
        
    def get_wavelengths(self,lenslet_num):
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
        npix = (distance * self.convert["mmtopx"]).astype(np.int) * self.density
        
        # Create a data array one hundred times as dense as the number of pixels
        #   This is the super dense array which will use the above interpolation
        dummy_lam = np.linspace(np.min(self.lams[use]),np.max(self.lams[use]),npix*100)
        
        # Interpolate along our really dense set of wavelengths to find all possible
        # illuminated pixel positions in this spectrum
        dummy_pts = np.array([fx(dummy_lam),fy(dummy_lam)])
        
        if self.density == 1:
            dummy_pts = dummy_pts.astype(np.int)
        else:
            dummy_pts = np.round(dummy_pts * self.density) / self.density
            
        # Remove any duplicate points. This does not do so in order, so we must
        # sort the array of unique points afterwards...
        unique_pts = np.array([np.array(x) for x in set(tuple(x) for x in dummy_pts.T)])
        
        # An array of distances to the origin of this spectrum, can be used to find wavelength
        # of light at each distance
        distance = np.array([np.sqrt(np.sum((x * self.convert["pxtomm"])-start)**2) for x in unique_pts])
        
        sorted_idx = np.argsort(distance)
        
        distance = distance[sorted_idx]
        points = unique_pts[sorted_idx]
        
        wl = self.spline(distance)
        
        return points,wl,np.diff(wl)
        
        
    def circle_kern(self,radius,size=0,normalize=False):
        """Generate a Circle Kernel"""
        # This allows us to set the size of the kernel arbitrarily
        if size < radius:
            size = int(radius)
        else:
            size = int(size)
        radius = int(radius)
        x, y = np.mgrid[-size:size+1, -size:size+1]
        d = np.sqrt(x**2.0 + y**2.0)
        v = (d <= radius).astype(np.float)
        if normalize:
            return v / np.sum(v)
        else:
            return v
        
    def gauss_kern(self,stdev,size=0,stdevy=None,sizey=0):
        """ Returns a normalized 2D gauss kernel array for convolutions """
        if size < (stdev**2.0):
            size = np.int(stdev**2.0)
        else:
            size = np.int(size)
        if not stdevy:
            stdevy = stdev
        if sizey < (stdevy**2.0):
            sizey = np.int(stdevy**2.0)
        else:
            sizey = np.int(sizey)
        
        x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
        g = np.exp(-(x**2.0/np.float(stdev**2.0)+y**2.0/np.float(stdevy**2.0)))
        
        return g / g.sum()
        
    
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
        self.linear = sp.interpolate.interp1d(WL,MM)
        # Generate the SPLINE function:
        self.spvars = sp.interpolate.splrep(MM,WL)
        self.spline = lambda x: sp.interpolate.splev(x,self.spvars)

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
        xs += (self.image_sz["mm"]/2)
        ys += (self.image_sz["mm"]/2)

        # Find the xs and ys that are not within 0.1 mm of the edge of the detector...
        ok = (xs > 0.1) & (xs < self.image_sz["mm"]-0.1) & (ys > 0.1) & (ys < self.image_sz["mm"]-0.1)
        ix, p1, p2, lams, xs, ys = ix[ok], p1[ok], p2[ok], lams[ok], xs[ok], ys[ok]
        # We remove these positions because we don't want to generate spectra for them.

        # This simply generates a list of all of the lenslets
        self.lenslets = np.unique(ix)
        
        # Convert the xs and ys to pixel positions
        xpix = np.round(xs * self.convert["mmtopx"],0).astype(np.int)
        ypix = np.round(ys * self.convert["mmtopx"],0).astype(np.int)
        
        # Determine the center of the whole system by finding the x position that is closest to 0,0 in pupil position
        cntix = np.argmin(p1**2 + p2**2)
        self.center = (xs[cntix] * self.convert["mmtopx"], ys[cntix] * self.convert["mmtopx"])
        
        self.ix, self.p1, self.p2, self.lams, self.xs, self.ys = ix, p1, p2, lams, xs, ys
        self.xpix, self.ypix = xpix, ypix
    
    def loadOpticsData(self,laspec,dispspec):
        """Loads an optical conversion based on the lenslet array spec and dispersion spec files provided"""
        
        self.LoadDispersionData(dispspec)
        self.LoadLensletData(laspec)
    
    def get_dense_image(self,lenslet,spectrum):
        """This function returns a dense image array, with flux placed into single pixels."""
        points, wl, deltawl = self.get_wavelengths(lenslet)
        # Some notes about this function:
        #  points:  These points are in normal CCD Pixel units, but are not integers,
        #           as they correspond to dense pixels on the overdense images
        #  wl:      A list of wavelengths that correspond to the points above
        #  deltawl: The np.diff(wl) result, which can be used to handle the flux per pixel
        
        # Create the spectra, by calling the spectra with wavelength
        # (in meters, note the conversion from wl, which was originally in microns)
        # Then use the deltawl to get the true amount of flux in each pixel
        radiance = spectrum(wl*1e-6) * self.gain
        flux = radiance[1,:-1] * deltawl
        
        # Take our points out. Note from the above that we multiply by the density in order to do this
        x,y = (points * self.density)[:-1].T.astype(np.int)
        
        # Get the way in which those points correspond to actual pixels. As such, this array of points should have duplicates
        xint,yint = points.T.astype(np.int)
        
        # Zero-adjust our x and y points. They will go into a fake subimage anyways, so we don't care
        # for now where they would be on the real image
        x -= np.min(x)
        y -= np.min(y)
        
        # Get the approximate size of our spectra
        xdist = np.max(x)-np.min(x)
        ydist = np.max(y)-np.min(y)
        
        # Convert this size into an integer number of pixels for our subimage. This makes it
        # *much* easier to register our sub-image to the master, larger pixel image
        xdist += (self.density - xdist % self.density)
        ydist += (self.density - ydist % self.density)
        
        # Move our x and y coordinates to the middle of our sub image by applying padding below each one.
        x += self.padding * self.density
        y += self.padding * self.density
        
        # Create our sub-image, using the x and y width of the spectrum, plus 2 padding widths.
        # Padding is specified in full-size pixels to ensure that the final image is an integer
        # number of full-size pixels across.
        img = np.zeros((xdist+2*self.padding*self.density,ydist+2*self.padding*self.density))
        
        # Place the spectrum into the sub-image
        img[x,y] = flux
        
        # Find the first (by teh flatten method) corner of the subimage, useful for placing the sub-image into the full image.
        corner = np.array([ xint[np.argmax(x)], yint[np.argmin(y)]]) - np.array([self.padding,self.padding])
        
        return img, corner
        
        
    def get_sub_image(self,lenslet,spectrum):
        """Returns a sub-image for a given lenslet"""
        # This function gets the dense sub-image with the spectrum placed in single dense pixels
        img,corner = self.get_dense_image(lenslet,spectrum)
        # Convolve this spectrum with an appropriate image of the telescope
        img2 = sp.signal.convolve(img,self.TELIMG,mode='same')
        # Convolve again with an appropriate PSF
        img3 = sp.signal.convolve(img2,self.PSFIMG,mode='same')
        # Bin the image back down to the final pixel size
        small = bin(img4,self.density).astype(np.int16)
        
        return small,corner,(img,img2,img3)
    
    
    def get_sub_image_fast(self,lenslet,spectrum):
        """Uses a single convolution, does not return intermediate steps"""
        # This function gets the dense sub-image with the spectrum placed in single dense pixels
        img,corner = self.get_dense_image(lenslet,spectrum)
        # Convolve with the PSF and Telescope Image simultaneously
        img2 = sp.signal.convolve(img,self.FINIMG,mode='same')
        # Bin the image back down to the final pixel size
        small = bin(img2,self.density).astype(np.int16)
        
        return small,corner
    
    def cache_sed_subimage(self,lenslet,spectrum):
        """Generates a sub image, and saves that result to this object. Should be thread-safe."""
        small, corner = self.get_sub_image_fast(lenslet,spectrum)
        label = "SUBIMG%d" % lenslet
        self.save(small,label)
        self.frame().metadata=dict(Lenslet=lenslet,Corner=corner)
    
    def place_cached_sed(self,lenslet,label,dlabel):
        """Places a cached SED Subimage"""
        slabel = "SUBIMG%d" % lenslet
        try:
            subframe = self.frame(slabel)
        except KeyError as e:
            raise SEDLimits(str(e))
        subimg = subframe()
        mlenslet = subframe.metadata['Lenslet']
        mcorner = subframe.metadata['Corner']
        if mlenslet != lenslet:
            raise ValueError("Lenslet Number Mismatch %d:%d for state %s in %s" % (lenslet,mlenslet,slabel,self))
        self.place(subimg,mcorner,label,dlabel)
    
    def place_sed(self,lenslet,spectrum,label,dlabel):
        """Place an SED based on a lenslet number and a spectrum"""
        small, corner = self.get_sub_image_fast(lenslet,spectrum)
        # Place uses the corner position to insert a sub-image into the master image
        self.place(small,corner,label,dlabel)
        
    def place(self,img,corner,label,dlabel):
        """Place the given AstroObject.AnalyticSpectra.AnalyticSpectrum onto the SEDMachine Image"""

        xstart = corner[0]
        xend = xstart + img.shape[0]
        ystart = corner[1]
        yend = ystart + img.shape[1]
        data = self.data(dlabel)
        
        if data.shape[0] < xend or data.shape[1] < yend:
            raise SEDLimits
            
        if xstart < 0 or ystart < 0:
            raise SEDLimits
            
        if xend < 0 or yend < 0:
            raise SEDLimits
            
        data[xstart:xend,ystart:yend] += img
        
        self.save(data,label)
        
    
    def generate_blank(self):
        """Generates a blank SEDMachine Image"""
        self.save(np.zeros((self.image_sz["px"],self.image_sz["px"])).astype(np.int16),"Blank")


    def crop(self,x,y,xsize,ysize=None):
        """Crops the provided image to twice the specified size, centered around the x and y coordinates provided."""
        if not ysize:
            ysize = xsize
        cropped = self.states[self.statename].data[x-xsize:x+xsize,y-ysize:y+ysize]
        LOG.debug("Cropped and Saved Image")
        self.save(cropped,"Cropped")
