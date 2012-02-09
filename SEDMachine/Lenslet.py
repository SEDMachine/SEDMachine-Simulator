# 
#  Lenslet.py
#  An object model of SED Machine Lenslets, containing most of the lenslet action logic.
#  SED
#  
#  Created by Alexander Rudy on 2012-01-31.
#  Copyright 2012 Alexander Rudy. All rights reserved.
#  Version 0.3
# 

import numpy as np
# import pyfits as pf
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.signal
import scipy.interpolate
import yaml

import shapely as sh
import shapely.geometry

import arpytools.progressbar

import os
import logging,logging.handlers
import time
import copy
import collections

import AstroObject
import AstroObject.AstroSimulator
from AstroObject.AstroCache import *
from AstroObject.AstroSpectra import SpectraObject
from AstroObject.AstroImage import ImageObject,ImageFrame
from AstroObject.AnalyticSpectra import BlackBodySpectrum, AnalyticSpectrum, FlatSpectrum
from AstroObject.Utilities import *


__version__ = getVersion(__name__)
__all__ = ["SEDLimits","Lenslet"]

class SEDLimits(Exception):
    """A Basic Error-Differentiation Class.
    This error is used to express the fact that the SEDModel has encountered a spectrum which can't be placed as some part of it falls outside of the limits of the SED system.
    """
    pass


class Lenslet(ImageObject):
    """An object-representation of a lenslet"""
    def __init__(self,xs,ys,xpixs,ypixs,p1s,p2s,ls,ix,config,caches):
        super(Lenslet, self).__init__()
        self.log = logging.getLogger("Instrument")
        self.num = ix
        self.xs = xs
        self.ys = ys
        self.points = np.array([xs,ys]).T
        self.xpixs = xpixs
        self.ypixs = ypixs
        self.pixs = np.array([xpixs,ypixs]).T
        self.ps = np.array([p1s,p2s]).T
        self.ls = np.array(ls)
        self.config = config
        self.Caches = caches        
        self.Caches.registerNPY("%dDispersion" % self.num,generate=self.return_dispersion,filename="Dispersion-%d.cache.npy" % self.num)
                
        self.dispersion = False
        self.checked = False
        self.passed = False
        self.traced = False
    
    def reset(self):
        """Reset Flag Variables"""
        self.dispersion = False
        self.checked = False
        self.passed = False
        self.traced = False
        
    
    def introspect(self):
        """Show all sorts of fun data about this lenslet"""
        STR  = "--Lenslet %(index)04d is %(valid)s\n" % {'index':self.num, 'valid': 'valid' if self.valid() else 'invalid'}
        STR += "|    x    |    y    |    xp    |    yp    |    p1    |    p2    |    wl    |\n"
        for xy,pixs,p,wl in zip(self.points,self.pixs,self.ps,self.ls):
            data = { 'x': xy[0], 'y': xy[1], 'pA': p[0], 'pB': p[1], 'wl': wl ,'pxA':pixs[0],'pxB':pixs[1]}
            STR += "|%(x) 9.6g|%(y) 9.6g|%(pxA) 10.6g|%(pxB) 10.6g|%(pA) 10.6g|%(pB) 10.6g|%(wl) 10.6g|\n" % data
        return STR
    
    def valid(self):
        """Returns true if this is a valid lenslet, false if it fails any of the tests"""
        if self.checked:
            return self.passed
        
        self.checked = True
        self.passed = False
        
        # Data consistency
        if len(self.points) != len(self.ps) or len(self.points) != len(self.ls) or len(self.points) != len(self.pixs):
            self.log.warning("Lenslet %d failed: The data had inconsistent points" % self.num)
            return self.passed
        
        # Data utility
        if len(self.points) < 3:
            self.log.debug("Lenslet %d failed: There were fewer than three data points" % self.num)
            return self.passed
        if np.any(self.pixs.flatten == 0):
            self.log.debug("Lenslet %d failed: Some (x,y) were exactly zero" % self.num)
            return self.passed
        
        # X distance calculation (all spectra should be roughly constant in x, as they are fairly well aligned)
        # NOTE: There really isn't a whole lot to this requriement
        dist = 30
        if np.any(np.abs(np.diff(self.xpixs)) > dist):
            self.log.debug("Lenslet %d failed: x distance was more than %d" % (self.num,dist))
            return self.passed
        
        # The spectrum should span some finite distance
        startix = np.argmin(self.ls)
        endix = np.argmax(self.ls)
        start = np.array([self.xs[startix],self.ys[startix]])
        end = np.array([self.xs[endix],self.ys[endix]])

        # Get the total length of the spectra
        self.distance = np.sqrt(np.sum(end-start)**2)
        
        if self.distance == 0:
            self.log.debug("Lenslet %d failed: The points have no separating distance" % self.num)
            return self.passed
        
        self.passed = True
        
        # Warnings about our data go here.
        if np.any(self.ls < 1e-12) or np.any(self.ls > 1e-3):
            self.log.warning("The wavelengths provided for lenslet %d appear as if they aren't SI units." % self.num)
            self.log.debug(npArrayInfo(self.ls,"Lenslet %d Wavelengths" % self.num))
                
        return self.passed
    
    def return_dispersion(self):
        """Return the dispersion"""
        self.find_dispersion()
        return self.dis
        
    def find_dispersion(self):
        """Find the dispersion (dense, pixel aligned wavelength values) for this lenslet"""
        assert self.valid(), "Lenslet must contain valid data."
        if self.dispersion:
            return self.dispersion
        
        # Interpolation to convert from wavelength to pixels.
        #   The accuracy of this interpolation is not important.
        #   Rather, it is used to find the pixels where the light will fall
        #   and is fed an array that is very dense, used on this dense interpolation
        #   and then binned back onto pixels. Thus it will be used to get a list
        #   of all illuminated pixels.
        fx = np.poly1d(np.polyfit(self.ls, self.xpixs, 2))
        fy = np.poly1d(np.polyfit(self.ls, self.ypixs, 2))
        
        # Find the starting and ending position of the spectra
        startix = np.argmin(self.ls)
        endix = np.argmax(self.ls)
        start = np.array([self.xs[startix],self.ys[startix]])
        end = np.array([self.xs[endix],self.ys[endix]])
        
        # Get the total length of the spectra
        distance = np.sqrt(np.sum(end-start)**2)
        
        # This should have been checked in the validity function.
        if distance == 0:
            raise SEDLimits
        
        # Find the length in units of (int) pixels
        npix = (distance * self.config["Instrument"]["convert"]["mmtopx"]).astype(np.int) * self.config["Instrument"]["density"]
        
        # Create a data array one hundred times as dense as the number of pixels
        #   This is the super dense array which will use the above interpolation
        superDense_lam = np.linspace(np.min(self.ls),np.max(self.ls),npix*100)
        
        # Interpolate along our really dense set of wavelengths to find all possible
        # illuminated pixel positions in this spectrum
        superDense_pts = np.array([fx(superDense_lam),fy(superDense_lam)])
        
        # Measure the distance along our really dense set of points
        superDense_interval = np.sqrt(np.sum(np.power(np.diff(superDense_pts,axis=1),2),axis=0))
        superDense_distance = np.cumsum(superDense_interval)
        
        # Adjust the density of our points. This rounds all values to only full pixel values.
        superDense_pts = np.round(superDense_pts * self.config["Instrument"]["density"]) / self.config["Instrument"]["density"]
        superDense_int = (superDense_pts * self.config["Instrument"]["density"]).astype(np.int)
        
        # We can identify unique points using the points when the integer position ratchets up or down.
        unique_x,unique_y = np.diff(superDense_int).astype(np.bool)
        
        # We want unique index to include points where either 'y' or 'x' ratchets up or down
        unique_idx = np.logical_or(unique_x,unique_y)
        
        # Remove any duplicate points. This does not do so in order, so we must
        # sort the array of unique points afterwards...
        unique_pts = superDense_pts[:,1:][:,unique_idx]
        
        # An array of distances to the origin of this spectrum, can be used to find wavelength
        # of light at each distance
        distance = superDense_distance[unique_idx] * self.config["Instrument"]["convert"]["pxtomm"]
        
        # Re sort everything by distnace along the trace.
        # Strictly, this shouldn't be necessary if all of the above functions preserved order.
        sorted_idx = np.argsort(distance)
        
        # Pull out sorted valuses
        distance = distance[sorted_idx]
        points = unique_pts[:,sorted_idx].T
        self.log.debug(npArrayInfo(points,"Points"))
        
        
        # Pull out the original wavelengths
        wl_orig = superDense_lam[unique_idx][sorted_idx]
        wl = wl_orig
        self.log.debug(npArrayInfo(wl,"Wavelengths"))
        # We are getting some odd behavior, where the dispersion function seems to not cover the whole
        # arc length and instead covers only part of it. This causes much of our arc to leave the desired
        # and available wavelength boundaries. As such, I'm disabling the more accurate dispersion mode.
        
        # Convert to wavelength space along the dispersion spline.
        # wl = self.spline(distance)
        xs,ys = points.T
        self.dxs = xs
        self.dys = ys
        self.dwl = wl
        self.drs = distance
        self.dis = np.array([xs,ys,wl,distance])
        self.dispersion = True
        
        self._debug_dispersion()
        
        return self.dispersion
                
    def get_trace(self,spectrum):
        """Returns a trace of this """
        
        if self.traced:
            return self.traced
        if not self.dispersion:
            self.dxs,self.dys,self.dwl,self.drs = self.Caches.get("%dDispersion" % self.num)
            self.dispersion = True
            
        # Variables taken from the dispersion calculation
        points = np.array([self.dxs,self.dys]).T
        deltawl = np.diff(self.dwl)
        
        
        # Take our points out. Note from the above that we multiply by the density in order to do this
        xorig,yorig = (points * self.config["Instrument"]["density"])[:-1].T.astype(np.int)
        x,y = (points * self.config["Instrument"]["density"])[:-1].T.astype(np.int)
        # Get the way in which those points correspond to actual pixels.
        # As such, this array of points should have duplicates
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
        xdist += (self.config["Instrument"]["density"] - xdist % self.config["Instrument"]["density"])
        ydist += (self.config["Instrument"]["density"] - ydist % self.config["Instrument"]["density"])
        
        # Move our x and y coordinates to the middle of our sub image by applying padding below each one.
        x += self.config["Instrument"]["padding"] * self.config["Instrument"]["density"]
        y += self.config["Instrument"]["padding"] * self.config["Instrument"]["density"]
        
        # Find the first (by the flatten method) corner of the subimage,
        # useful for placing the sub-image into the full image.
        corner = np.array([ xint[np.argmax(x)], yint[np.argmin(y)]])
        self.log.debug("Corner Position in Integer Space: %s" % corner)
        corner *= self.config["Instrument"]["density"]
        realcorner = np.array([ xorig[np.argmax(x)], yorig[np.argmin(y)]])
        offset = corner - realcorner
        corner /= self.config["Instrument"]["density"]
        self.log.debug("Corner Position Offset in Dense Space: %s" % (offset))
        if self.log.getEffectiveLevel() <= logging.DEBUG:
            with open(self.config["Dirs"]["Partials"]+"Instrument-Offsets.dat",'a') as handle:
                np.savetxt(handle,offset)
        corner -= np.array([-self.config["Instrument"]["padding"],self.config["Instrument"]["padding"]])
        
        x += offset[0]
        y += offset[1]
        
        # Create our sub-image, using the x and y width of the spectrum, plus 2 padding widths.
        # Padding is specified in full-size pixels to ensure that the final image is an integer
        # number of full-size pixels across.
        xsize = xdist+2*self.config["Instrument"]["padding"]*self.config["Instrument"]["density"]
        ysize = ydist+2*self.config["Instrument"]["padding"]*self.config["Instrument"]["density"]
        
        # Calculate the resolution inherent to the pixels asked for
        WLS = self.dwl
        DWL = np.diff(WLS) 
        WLS = WLS[:-1]
        RS = WLS/DWL
        
        self.log.debug("Scailing by %g" % self.config["Instrument"]["gain"])
        radiance = spectrum(wavelengths=WLS,resolution=RS) 
        radiance *= self.config["Instrument"]["gain"]
        self.log.debug(npArrayInfo(radiance,"Generated Spectrum"))
        self.log.debug(npArrayInfo(deltawl,"DeltaWL Rescaling"))
        flux = radiance[1,:] * deltawl
        self.log.debug(npArrayInfo(flux,"Final Flux"))
         
        if self.log.getEffectiveLevel() <= logging.DEBUG:
            np.savetxt(self.config["Dirs"]["Partials"]+"Instrument-Subimage-Values.dat",np.array([x,y,self.dwl[:-1],deltawl,flux]).T)
        
        self.txs = x
        self.tys = y
        self.tfl = flux
        self.twl = self.dwl[:-1]
        self.subshape = (xsize,ysize)
        self.subcorner = corner
        self.traced = True
        
        return self.traced
        
    def place_trace(self,get_psf):
        """Place the trace on the image"""
        
        img = np.zeros(self.subshape)
        
        for x,y,wl,flux in zip(self.txs,self.tys,self.twl,self.tfl):
            psf = get_psf(wl)
            tiny_image = np.zeros(psf.shape) + psf * flux
            tl_corner = [ x - tiny_image.shape[0]/2.0, y - tiny_image.shape[0]/2.0 ]
            br_corner = [ x + tiny_image.shape[0]/2.0, y + tiny_image.shape[0]/2.0 ]
            self.log.debug("Corner of tiny image is %s" % (tl_corner))
            img[tl_corner[0]:br_corner[0],tl_corner[1]:br_corner[1]] += tiny_image
        self.log.debug(npArrayInfo(img,"DenseSubImage"))
        self.save(img,"Raw Spectrum")
        
        
    def _debug_dispersion(self):
        """Debugging for the dispersion process"""
        if self.log.getEffectiveLevel() > logging.DEBUG or not self.config["Plot"]:
            return
        self.log.debug("Plotting Wavelength Information")
        # This graph shows the change in distance along arc per pixel.
        # The graph should produce all data points close to each other, except a variety of much lower
        # data points which are caused by the arc crossing between pixels.
        plt.figure()
        plt.title("$\Delta$Distance Along Arc")
        plt.xlabel("x (px)")
        plt.ylabel("$\Delta$Distance along arc (px)")
        plt.plot(self.dys[:-1],np.diff(distance) * self.config["Instrument"]["convert"]["mmtopx"],'g.')
        plt.savefig("%sInstrument-%04d-Delta-Distances%s" % (self.config["Dirs"]["Partials"],self.num,self.config["plot_format"]))
        plt.clf()
    
    def rotate(self,point,angle,origin=None):
        """Rotate a given point by the provided angle around the origin given"""
        if origin == None:
            origin = np.array([0,0])
        pA = np.matrix((point - origin))
        R = np.matrix([[np.cos(angle),-1*np.sin(angle)],[np.sin(angle),np.cos(angle)]]).T
        pB = pA * R
        return np.array(pB + origin)[0]
    
    def make_hexagon(self):
        """Generates a hexagonal polygon for this lenslet, to be used for area calculations"""
        radius = 0.25e-2
        angle = np.pi/3.0
        rotation = 20.0 * (np.pi/180.0) #Absolute angle of rotation for hexagons
        A = self.rotate(self.ps + np.array([radius,0]),rotation,self.ps)
        points = [A]
        for i in range(6):
            A = self.rotate(A,angle,self.ps)
            points.append(A)
        self.shape = sh.geometry.Polygon(tuple(points))

        
        

