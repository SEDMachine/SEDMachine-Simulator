# 
#  Objects.py
#  An object model of SED Machine Lenslets, containing most of the lenslet action logic.
#  SED
#  
#  Created by Alexander Rudy on 2012-01-31.
#  Copyright 2012 Alexander Rudy. All rights reserved.
#  Version 0.3.9-p4
# 

import numpy as np
import pyfits as pf
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm # Python Color Systems
from matplotlib.colors import Normalize # Math Normailzation


import scipy.signal
import scipy.interpolate
import yaml

import shapely as sh
import shapely.geometry

import os
import logging,logging.handlers
import time
import copy
import collections
import gc

import AstroObject
import AstroObject.AstroSimulator
from AstroObject.AstroCache import *
from AstroObject.AstroSpectra import SpectraStack
from AstroObject.AstroImage import ImageStack,ImageFrame
from AstroObject.AnalyticSpectra import BlackBodySpectrum, AnalyticSpectrum, FlatSpectrum, InterpolatedSpectrum
from AstroObject.util import getVersion, npArrayInfo


__version__ = getVersion()
__all__ = ["SEDLimits","Lenslet","SubImage","SourcePixel"]

class SEDLimits(Exception):
    """A Basic Error-Differentiation Class.
    This error is used to express the fact that the SEDModel has encountered a spectrum which can't be placed as some part of it falls outside of the limits of the SED system.
    """
    pass

class SubImage(ImageFrame):
    """A custom frame type for SEDMachine Sub-Images"""
    def __init__(self, array, label, header=None, metadata=None):
        super(SubImage, self).__init__(array,label,header,metadata)
        self.lensletNumber = 0
        self.corner = [0,0]
        self.configHash = hash(0)
        self.spectrum = "NO SPEC"

    log = logging.getLogger("SEDMachine")
    
        
    def sync_header(self):
        """Synchronizes the header dictionary with the HDU header"""
        # assert self.label == self.header['SEDlabel'], "Label does not match value specified in header: %s and %s" % (self.label,self.header['SEDlabel'])
        self.configHash = self.header['SEDconf']
        self.corner = [self.header['SEDcrx'],self.header['SEDcry']]
        self.spectrum = self.header['SEDspec']
        self.lensletNumber = self.header['SEDlens']
    
    def __hdu__(self,primary=False):
        """Retruns an HDU which represents this frame. HDUs are either ``pyfits.PrimaryHDU`` or ``pyfits.ImageHDU`` depending on the *primary* keyword."""
        if primary:
            self.log.log(5,"Generating a primary HDU for %s" % self)
            HDU = pf.PrimaryHDU(self())
        else:
            self.log.log(5,"Generating an image HDU for %s" % self)
            HDU = pf.ImageHDU(self())
        return HDU
        
    def __setheader__(self,HDU):
        """Sets the header object"""
        HDU.header.update('object',self.label)
        HDU.header.update('SEDlabel',self.label)
        HDU.header.update('SEDconf',self.configHash)
        HDU.header.update('SEDcrx',self.corner[0])
        HDU.header.update('SEDcry',self.corner[1])
        HDU.header.update('SEDspec',self.spectrum)
        HDU.header.update('SEDlens',self.lensletNumber)
        return HDU        
    
    def __show__(self):
        """Plots the image in this frame using matplotlib's ``imshow`` function. The color map is set to an inverted binary, as is often useful when looking at astronomical images. The figure object is returned, and can be manipulated further.
        
        .. Note::
            This function serves as a quick view of the current state of the frame. It is not intended for robust plotting support, as that can be easily accomplished using ``matplotlib``. Rather, it attempts to do the minimum possible to create an acceptable image for immediate inspection.
        """
        self.log.debug("Plotting %s using matplotlib.pyplot.imshow" % self)
        figure = plt.imshow(self())
        figure.set_cmap('binary_r')
        return figure
    
    @classmethod
    def __read__(cls,HDU,label):
        """Attempts to convert a given HDU into an object of type :class:`ImageFrame`. This method is similar to the :meth:`__save__` method, but instead of taking data as input, it takes a full HDU. The use of a full HDU allows this method to check for the correct type of HDU, and to gather header information from the HDU. When reading data from a FITS file, this is the prefered method to initialize a new frame.
        """
        cls.log.debug("Attempting to read as %s" % cls)
        if not isinstance(HDU,(pf.ImageHDU,pf.PrimaryHDU)):
            msg = "Must save a PrimaryHDU or ImageHDU to a %s, found %s" % (cls.__name__,type(HDU))
            raise AbstractError(msg)
        if not isinstance(HDU.data,np.ndarray):
            msg = "HDU Data must be %s for %s, found data of %s" % (np.ndarray,cls.__name__,type(HDU.data))
            raise AbstractError(msg)
        try:
            Object = cls(HDU.data,label,header=HDU.header)
        except AssertionError as AE:
            msg = "%s data did not validate: %s" % (cls.__name__,AE)
            raise AbstractError(msg)
        cls.log.debug("Created %s" % Object)
        Object.sync_header()
        return Object


class Lenslet(ImageStack):
    """An object-representation of a lenslet. Takes approximately all of the data we know about each lenslet.
    
    :param xs: Array of camera-center x positions
    :param ys: Array of camera-center y positions
    :param p1s: Array of pupil x positions
    :param p2s: Array of pupil y positions
    :param ls: Array of wavelengths
    :param ix: Index of this lenslet
    :param config: Configuration object
    :param caches: Cache object
    
    """
    def __init__(self,p1s,p2s,ls,ix,xcs, ycs, xls, yls, xas, yas, xbs, ybs, rs,config,caches):
        super(Lenslet, self).__init__()
        self.dataClasses = [SubImage]
        self.log = logging.getLogger("SEDMachine")
        self.config = config
        self.num = ix
        
        # Center x and y positions
        self.xcs =  xcs + (self.config["Instrument"]["image"]["size"]["mm"]/2)
        self.ycs =  ycs + (self.config["Instrument"]["image"]["size"]["mm"]/2)
        self.xas =  xas + (self.config["Instrument"]["image"]["size"]["mm"]/2)
        self.yas =  yas + (self.config["Instrument"]["image"]["size"]["mm"]/2)
        self.xbs =  xbs + (self.config["Instrument"]["image"]["size"]["mm"]/2)
        self.ybs =  ybs + (self.config["Instrument"]["image"]["size"]["mm"]/2)
        self.xls =  xls + (self.config["Instrument"]["image"]["size"]["mm"]/2)
        self.yls =  yls + (self.config["Instrument"]["image"]["size"]["mm"]/2)
        
        
        
        self.points = np.array([self.xcs,self.ycs]).T
        
        # Convert the xs and ys to pixel positions
        self.xpixs = np.round(self.xcs * self.config["Instrument"]["convert"]["mmtopx"],0).astype(np.int)
        self.ypixs = np.round(self.ycs * self.config["Instrument"]["convert"]["mmtopx"],0).astype(np.int)
        self.pixs = np.array([self.xpixs,self.ypixs]).T
        self.ps = np.array([p1s,p2s]).T
        self.ls = np.array(ls)
                
        self.dispersion = False
        self.checked = False
        self.passed = False
        self.traced = False
        self.spectrum = FlatSpectrum(0.0)
        
    
    def reset(self):
        """Reset Flag Variables for this lenslet. Deletes any calculated varaibles."""
        self.checked = False
        self.passed = False
        
        if self.dispersion:
            self.dispersion = False
            del self.dxs
            del self.dys
            del self.dwl
            del self.drs
            del self.dis
        
        if self.traced:
            self.traced = False
            del self.txs
            del self.tys
            del self.trd
            del self.tfl
            del self.twl
            del self.tdw
            del self.trs
            del self.subshape
            del self.subcorner
            
        gc.collect()
        
        
    
    def introspect(self):
        """Show all sorts of fun data about this lenslet.
        
        ::
            
            >>> lenslet.introspect()
            --Lenslet 0023 is valid
            |    x    |    y    |    xp    |    yp    |    p1    |    p2    |    wl    |
            |  37.5052|  39.0077|      2778|      2889|   -0.0692|   -0.1593|   3.7e-07|
            |  37.8836|  40.6228|      2806|      3009|   -0.0692|   -0.1593|  6.45e-07|
            |  37.9199|  41.4872|      2809|      3073|   -0.0692|   -0.1593|   9.2e-07|
            
        
        """
        STR  = "--Lenslet %(index)04d is %(valid)s\n" % {'index':self.num, 'valid': 'valid' if self.valid() else 'invalid'}
        STR += "|    x    |    y    |    xp    |    yp    |    p1    |    p2    |    wl    |\n"
        for xy,pixs,p,wl in zip(self.points,self.pixs,self.ps,self.ls):
            data = { 'x': xy[0], 'y': xy[1], 'pA': p[0], 'pB': p[1], 'wl': wl ,'pxA':pixs[0],'pxB':pixs[1]}
            STR += "|%(x) 9.6g|%(y) 9.6g|%(pxA) 10.6g|%(pxB) 10.6g|%(pA) 10.6g|%(pB) 10.6g|%(wl) 10.6g|\n" % data
        return STR
    
    def valid(self,strict=True):
        """Returns true if this is a valid lenslet, false if it fails any of the tests.
        
        :param strict: In not strict mode, validator will let many lenslets which are not well-formed through the system.
        :returns: bool
        
        Checks performed:
        
        - Wavelengths and Points all have the same number of entries.
        - We have at least three entries.
        - No entry's pixel position is exactly 0.
        - The dispersion distance along the x-axis is less than 30 pixels.
        - The distance between the start and end of the spectrum is non-zero.
        - The lenslet positions are not within some configured tolerance of the edge of the spectrum.
        - Warning about the units for wavelengths.
        
        
        
        """
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
            if strict:
                return self.passed
        if np.any(self.pixs.flatten == 0):
            self.log.debug("Lenslet %d failed: Some (x,y) were exactly zero" % self.num)
            if strict:
                return self.passed
        
        # X distance calculation (all spectra should be roughly constant in x, as they are fairly well aligned)
        # NOTE: There really isn't a whole lot to this requriement
        dist = 30
        if np.any(np.abs(np.diff(self.xpixs)) > dist):
            self.log.debug("Lenslet %d failed: x distance was more than %d" % (self.num,dist))
            if strict:
                return self.passed
        
        # The spectrum should span some finite distance
        startix = np.argmin(self.ls)
        endix = np.argmax(self.ls)
        start = np.array([self.xcs[startix],self.ycs[startix]])
        end = np.array([self.xcs[endix],self.ycs[endix]])

        # Get the total length of the spectra
        self.distance = np.sqrt(np.sum(end-start)**2)
        
        if self.distance == 0:
            self.log.debug("Lenslet %d failed: The points have no separating distance" % self.num)
            return self.passed
        
        # Find the xs and ys that are not within 0.1 mm of the edge of the detector...
        padding = self.config["Instrument"]["image"]["pad"]["mm"]
        if not ((self.xcs > 0.1) & (self.xcs < self.config["Instrument"]["image"]["size"]["mm"]-padding) & (self.ycs > padding) & (self.ycs < self.config["Instrument"]["image"]["size"]["mm"]-padding)).any():
            self.log.debug("Lenslet %d failed: The points are too close to the image edge" % self.num)
            if strict:
                return self.passed
        
        self.passed = True
        
        # Warnings about our data go here.
        if np.any(self.ls < 1e-12) or np.any(self.ls > 1e-3):
            self.log.warning("The wavelengths provided for lenslet %d appear as if they aren't SI units." % self.num)
            self.log.debug(npArrayInfo(self.ls,"Lenslet %d Wavelengths" % self.num))
                
        return self.passed
    
    def return_dispersion(self):
        """Return the dispersion values.
        
        :returns: Array of [xcs,ycs,wl,distance]
        
        See :meth:`find_dispersion` for the calculation performed.
        """
        self.find_dispersion()
        return self.dis
        
    def find_dispersion(self):
        """Find the dispersion (dense, pixel aligned wavelength values) for this lenslet.
        
        To calculate dispersion, we first create an interpolation from (wavelength) -> (xpix) and (wavelength) -> (ypix). Then, using a large array of possible wavelengths (100x the number of oversampled pixel positions) we create arrays of possible overdense x and y pixel positions. These arrays are then truncated to contain only integer (overdense) pixel positions. We then calculate the arc-distance to each of these pixel positions from the start (lowest wavelength) of the spectrum. Taking only unique points along the arc-distance, we find a list of all of the unique x and y pixel positions which are illuminated by the spectrum in the over-dense sample space. This array, along with thier corresponding wavelengths and arc-distances, are stored for later use.
        
        **Variables which are used**:
        
        :var xcs: x-camera positions (center of image) in mm
        :var ycs: y-camera positions (center of image) in mm
        :var ls: wavelengths
        :var xpixs: x-camera positions (center of image) in px (integer)
        :var ypixs: y-camera positions (center of image) in px (integer)
        
        
        **Variables which are set**:
        
        :var dxs: x-camera-positions of each illuminated oversampled pixel in px (integer o-px)
        :var dys: y-camera-positions of each illuminated oversampled pixel in px (integer o-px)
        :var dwl: wavelength of each illuminated oversampled pixel in meters
        :var drs: arc-distance along spectrum in mm
        :var dis: array of ``[dxs,dys,dwl,drs]``
        :var dispersion: boolean True
        
        """
        assert self.valid(), "Lenslet must contain valid data."
        if self.dispersion:
            return self.dispersion
        
        if self.config["Instrument"]["Tel"]["ellipse"]:
            # Find ellipse major and minor axis from given data.
            self.a = np.sqrt((self.xcs - self.xas)**2.0 + (self.ycs-self.yas)**2.0) * self.config["Instrument"]["convert"]["mmtopx"] * self.config["Instrument"]["density"]
            self.b = np.sqrt((self.xcs - self.xbs)**2.0 + (self.ycs-self.ybs)**2.0) * self.config["Instrument"]["convert"]["mmtopx"] * self.config["Instrument"]["density"]
            top = self.xcs - self.xas
            bot = self.ycs - self.yas
            bot[np.logical_and(top == 0,bot == 0)] = 1.0
            self.alpha = np.arctan(top/bot)
            self.fa = np.poly1d(np.polyfit(self.ls, self.a, self.config["Instrument"]["Tel"]["dispfitorder"]))
            self.fb = np.poly1d(np.polyfit(self.ls, self.b, self.config["Instrument"]["Tel"]["dispfitorder"]))
            self.falpha = np.poly1d(np.polyfit(self.ls, self.alpha, self.config["Instrument"]["Tel"]["dispfitorder"]))
            

            
        
        # Interpolation to convert from wavelength to pixels.
        #   The accuracy of this interpolation is not important.
        #   Rather, it is used to find the pixels where the light will fall
        #   and is fed an array that is very dense, used on this dense interpolation
        #   and then binned back onto pixels. Thus it will be used to get a list
        #   of all illuminated pixels.
        fx = np.poly1d(np.polyfit(self.ls, self.xpixs, self.config["Instrument"]["dispfitorder"]))
        fy = np.poly1d(np.polyfit(self.ls, self.ypixs, self.config["Instrument"]["dispfitorder"]))
        
        # Find the starting and ending position of the spectra
        startix = np.argmin(self.ls)
        endix = np.argmax(self.ls)
        start = np.array([self.xcs[startix],self.ycs[startix]])
        end = np.array([self.xcs[endix],self.ycs[endix]])
        
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
        
        return self.dispersion
                
    def get_trace(self,spectrum):
        """Returns a trace of this spectrum. The trace will contain x and y over-dense pixel positions, flux values for each of those illuminated pixels, instantaneous resolution at each pixel, and wavelength of each pixel. The trace also determines the corners of the spectrum, and saves those corner positions with ample padding in the image.
        
        To perform the trace, we first get the oversampled point positions in o-px from the ``dxs`` and ``dys`` variables. These points are saved as both ``x,y`` and ``xorig,yorig``, for later use. ``xorig,yorig`` are stored to save an unmodified copy of the points. The ``xint,yint`` variables are used to store the integer (in px) positions of each ``x,y`` pair, essentially, thier containing camera pixel. ``x,y`` are then zeroed, such that (0,0) is the upper-right corner of the spectrum to be inserted. We then calculate the size of the subimage (in ``xdist,ydist``) and adjust this size so that it is an integer number of camera pixels (px) across. This makes the binning much easier later. The pixels ``x,y`` are then padded to provide space for the PSF to be applied on all sides of the single-pixel spectrum. 
        
        We then find the corner of the of the subimage. First, the corner will generally be extracted as the position with a minimum in the ``y`` direction and a maximum in the ``x`` direction. We first find the corner's position in integer camrea-px space (i.e. from ``xint,yint``, stored as ``corner``), and the corner's position in integer o-px space (i.e. from ``xorig,yorig``, stored as ``realcorner``). Converting the camera's position in integer camera-px space, we take the difference between the two corners as the ``offset``. This is the shift we must insert into ``x,y`` in order to ensure that the corner of our subimage will line up with the corner of a binned pixel. Next we add padding distances into the ``corner`` position. Finally, we use the offset to move the ``x,y`` positions to account for aligning the corner of the subimage with the corner of a full camera pixel. I will make a diagram to explain all of this shortly. Next, we add the padding values into ``xdist,ydist`` to get the full size of the subimage in o-px.
        
        We are now ready to extract flux values from our spectrum. This is done using the ``dwl`` values as the wavelength values to sample at. Using the ``dwl`` values, we also calculate an effective sampling resolution, which is used to resample the spectra. Feeding both of these, we compute the flux of the spectrum at each pixel position. This computation is not described in this function, but in a separate location in the documentation.
        
        After computing the flux at each pixel, we save the ``x,y`` indicies of each pixel, the flux for that pixel, the wavelength of that pixel, and the shape and corner of the subimage for later use.
        
        
        **Variables which are used**:
        
        :var dxs: x-camera-positions of each illuminated oversampled pixel in px
        :var dys: y-camera-positions of each illuminated oversampled pixel in px
        :var dwl: wavelength of each illuminated oversampled pixel in meters
        
        **Variables which are set**:
        
        :var txs: x-subimage-indicies of each illuminated oversampled pixel in o-px
        :var tys: y-subimage-indicies of each illuminated oversampled pixel in o-px
        :var tfl: flux of each pixel in counts
        :var twl: wavelength of each pixel in meters
        :var tdw: delta wavelength covered by each pixel in meters
        :var trs: sampling resolution of each pixel
        :var subshape: shape of subimage to contain spectrum
        :var subcorner: corner of subimage to contain spectrum in px
        :var spectrum: the spectrum object used for flux
        :var traced: bool True
        
        """
        
        if self.traced:
            return self.traced
            
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
            with open("%(Partials)s/Instrument-Offsets.dat" % self.config["Dirs"],'a') as handle:
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
        
        
        # Call and evaluate the spectrum
        self.log.debug(npArrayInfo(WLS,"Calling Wavelength"))
        self.log.debug(npArrayInfo(RS,"Calling Resolution"))

        wl,flux = spectrum(wavelengths=WLS,resolution=RS) 
                
        self.log.debug(npArrayInfo(flux,"Final Flux"))
        self.log.debug(npArrayInfo(RS,"Saving Resolution"))
        
        self.txs = x
        self.tys = y
        self.tfl = flux
        self.twl = WLS
        self.tdw = DWL
        self.trs = RS
        self.subshape = (xsize,ysize)
        self.subcorner = corner
        self.spectrum = spectrum
        self.traced = True
        
        return self.traced
        
    def place_trace(self,get_conv):
        """Place the trace on the subimage.
        
        First a blank image is created using the ``subshape`` variable as a template. Then, iterating through each point, we get the convolved PSF and telescope image for that point (called the "convolution"). The convolution can vary by wavelenght, and can have other variables which are sytem dependent. The convolution is multiplied by the flux value for that point. The corner of the convolution is then calculated (the top-left and bottom right corners are actually calculated) so that the image can be insterted as a 'flattened' array. Each image is inserted by addition, adding on to the image already in place.
        
        **Variables which are used**:
        
        :var txs: x-subimage-indicies of each illuminated oversampled pixel in o-px
        :var tys: y-subimage-indicies of each illuminated oversampled pixel in o-px
        :var tfl: flux of each pixel in counts
        :var twl: wavelength of each pixel in meters
        :var subshape: shape of subimage to contain spectrum
        :var subcorner: corner of subimage to contain spectrum in px
        
        **Variables which are set**:
        
        :var frame: A frame labeled "Raw Spectrum" with the oversampled spectrum.
        
        """
        
        img = np.zeros(self.subshape)
        
        for x,y,wl,flux in zip(self.txs,self.tys,self.twl,self.tfl):
            if self.config["Instrument"]["Tel"]["ellipse"]:
                a = self.fa(wl)
                b = self.fb(wl)
                rot = self.falpha(wl)
                conv = get_conv(wl,a,b,rot)
            else:
                conv = get_conv(wl)
            tiny_image = conv * flux
            tl_corner = [ x - tiny_image.shape[0]/2.0, y - tiny_image.shape[0]/2.0 ]
            br_corner = [ x + tiny_image.shape[0]/2.0, y + tiny_image.shape[0]/2.0 ]
            img[tl_corner[0]:br_corner[0],tl_corner[1]:br_corner[1]] += tiny_image
        self.log.debug(npArrayInfo(img,"DenseSubImage"))
        self["Raw Spectrum"] = img
        frame = self.frame()
        frame.lensletNumber = self.num
        frame.corner = self.subcorner
        frame.configHash = hash(str(self.config.extract()))
    
    def write_subimage(self):
        """Writes a subimage to file and then clears the subimage from this lenslet's memory. The subimage is written to a file with a name formatted as ``%(Caches)s/Subimage-%(num)4d.fits`` where the fields to the format string are relatively self explanatory.
        
        Subimage writes will clobber old images.
        """
        fileName = "%(Caches)s/Subimage-%(num)04d%(ext)s" % dict(num=self.num,ext=".fits",**self.config["Dirs"])
        if os.access(fileName,os.F_OK):
            os.remove(fileName)
        self.write(fileName,primaryFrame="Raw Spectrum",clobber=True)
        self.clear()
        
    def read_subimage(self):
        """Read a subimage from file and set the lenslet number and corner from data contained in the file. The filenames for subimages conform to the format ``%(Caches)s/Subimage-%(num)4d.fits``, similar to :meth:`write_subimage`. Frames are :class:`SubImage` classes so that they can safely store and retrieve their lenslet numbers and corner positions."""
        self.read("%(Caches)s/Subimage-%(num)04d%(ext)s" % dict(num=self.num,ext=".fits",**self.config["Dirs"]))
        frame = self.frame()
        self.num = frame.lensletNumber
        self.subcorner = frame.corner
        
    def bin_subimage(self):
        """Bin the selected subimage using the :meth:`bin` function, and binning based on the configured density. This function also sets the final data type as ``np.int16``.
        
        """
        self["Binned Spectrum"] = self.bin(self.data(),self.config["Instrument"]["density"])
    
    def plot_raw_data(self):
        """Save a plot figure for raw-data from the lenslet.
        
        Creates 1 plot of y-position vs. wavelegnth for the given data points for each lenslet. Plot is saved as ``%(Partials)s/Lenslet-%(num)04d-WL%(ext)s``.
        
        **Variables which are used**:
        
        :var ls: Wavelengths in meters
        :var ycs: y-camera positions (center of image) in mm
        
        """
        plt.clf()
        plt.plot(self.ls*1e6,self.ycs,".",linestyle='-')
        plt.title("$\lambda$ along y-axis (%d)" % self.num)
        plt.xlabel("Wavelength ($\mu m$)")
        plt.ylabel("Y-position ($mm$)")
        plt.savefig("%(Partials)s/Lenslet-%(num)04d-WL%(ext)s" % dict(num=self.num, ext=self.config["Plots"]["format"],**self.config["Dirs"]))
        plt.clf()
        plt.plot(self.ls*1e6,self.ycs,".",linestyle='-',label="Center")
        plt.plot(self.ls*1e6,self.yas,".",linestyle='-',label="Major Axis")
        plt.title("$\lambda$ along y-axis (%d)" % self.num)
        plt.xlabel("Wavelength ($\mu m$)")
        plt.ylabel("Y-position ($mm$)")
        plt.legend()
        plt.savefig("%(Partials)s/Lenslet-%(num)04d-WL-ya%(ext)s" % dict(num=self.num, ext=self.config["Plots"]["format"],**self.config["Dirs"]))
        plt.clf()
        plt.subplot(211)
        plt.plot(self.ls*1e6,np.abs(self.xcs-self.xbs)*self.config["Instrument"]["convert"]["mmtopx"],".",linestyle='-',label="XB")
        plt.plot(self.ls*1e6,np.abs(self.ycs-self.yas)*self.config["Instrument"]["convert"]["mmtopx"],".",linestyle='-',label="YA")
        plt.title("$\lambda$ along y-axis (%d)" % self.num)
        plt.xlabel("Wavelength ($\mu m$)")
        plt.ylabel("$\Delta$Y-position ($px$)")
        plt.legend()
        plt.axis(expandLim(plt.axis()))

        plt.subplot(212)
        plt.plot(self.ls*1e6,np.abs(self.xcs-self.xas)*self.config["Instrument"]["convert"]["mmtopx"],".",linestyle='-',label="XA")
        plt.plot(self.ls*1e6,np.abs(self.ycs-self.ybs)*self.config["Instrument"]["convert"]["mmtopx"],".",linestyle='-',label="YB")
        plt.axis(expandLim(plt.axis()))
        plt.xlabel("Wavelength ($\mu m$)")
        plt.ylabel("$\Delta$Y-position ($px$)")
        plt.legend()
        plt.savefig("%(Partials)s/Lenslet-%(num)04d-WL-dy%(ext)s" % dict(num=self.num, ext=self.config["Plots"]["format"],**self.config["Dirs"]))
        plt.clf()
        plt.plot(self.ls*1e6,self.fa(self.ls),".",linestyle='-')
        plt.title("$\lambda$ along y-axis (%d)" % self.num)
        plt.xlabel("Wavelength ($\mu m$)")
        plt.ylabel("Major-Axis ($px$)")
        plt.savefig("%(Partials)s/Lenslet-%(num)04d-WL-a%(ext)s" % dict(num=self.num, ext=self.config["Plots"]["format"],**self.config["Dirs"]))
        plt.clf()
        
    def plot_ellipses(self):
        """docstring for plot_ellipses"""
        plt.plot(self.ls*1e6,np.abs(self.ycs-self.yas)*self.config["Instrument"]["convert"]["mmtopx"],".",linestyle='-',label="YA")


    def plot_rotation(self):
        """docstring for plot_ellipses"""
        plt.plot(self.ls*1e6,self.falpha(self.ls) * 180.0 / np.pi,".",linestyle='-',label="YA")
                
        
    def plot_dispersion(self):
        """Two plots which show results from dispersion calculations.
        
        1. Plot of arc-distance (delta) in px vs. x in pixels. ``%(Partials)s/Instrument-%(num)04d-Delta-Distances%(ext)s``
        
        2. Plot of Wavelength vs. effective sampling resolution. ``%(Partials)s/Instrument-%(num)04d-WL%(ext)s``
        
        **Variables which are used**:
        
        :var dys: y-camera-positions of each illuminated oversampled pixel in px (integer o-px)
        :var dwl: wavelength of each illuminated oversampled pixel in meters
        :var drs: arc-distance along spectrum in mm
        
        """
        assert self.dispersion
        # This graph shows the change in distance along arc per pixel.
        # The graph should produce all data points close to each other, except a variety of much lower
        # data points which are caused by the arc crossing between pixels.
        plt.figure()
        plt.title("$\Delta$Distance Along Arc (%d)" % self.num)
        plt.xlabel("x (px)")
        plt.ylabel("$\Delta$Distance along arc (px)")
        plt.plot(self.dys[:-1],np.diff(self.drs) * self.config["Instrument"]["convert"]["mmtopx"],'g.')
        plt.savefig("%(Partials)s/Instrument-%(num)04d-Delta-Distances%(ext)s" % dict(num=self.num, ext=self.config["Plots"]["format"],**self.config["Dirs"]))
        plt.clf()
        plt.plot(self.drs,self.dwl*1e6,"g.")
        plt.title("$\lambda$ for each pixel (%d)" % self.num)
        plt.ylabel("Wavelength ($\mu m$)")
        plt.xlabel("Resolution $R = \\frac{\lambda}{\Delta\lambda}$ per pixel")
        plt.savefig("%(Partials)s/Instrument-%(num)04d-WL%(ext)s" % dict(num=self.num, ext=self.config["Plots"]["format"],**self.config["Dirs"]))
        plt.clf()
        
        
        
    def plot_spectrum(self):
        """One plot showing the flux of the spectrum for this lenslet.
        
        Plot shows Flux (in counts) vs. Wavelength. ``%(Partials)s/Instrument-%(num)04d-Flux%(ext)s``
        
        **Variables which are used**:
        
        :var tfl: flux of each pixel in counts
        :var twl: wavelength of each pixel in meters
        
        """
        assert self.traced
        self.log.debug(npArrayInfo(self.twl*1e6,"Wavlength"))
        self.log.debug(npArrayInfo(self.tfl,"Flux"))
        plt.clf()
        plt.semilogy(self.twl*1e6,self.tfl,"b.")
        xmin,xmax,ymin,ymax = plt.axis()
        if ymin < 1e-4:
            ymin = 1e-4
        plt.axis((xmin,xmax,ymin,ymax))
        plt.title("Generated, Fluxed Spectra (%d)" % self.num)
        plt.xlabel("Wavelength ($\mu m$)")
        plt.ylabel("Flux (Electrons)")
        plt.savefig("%(Partials)s/Instrument-%(num)04d-Flux%(ext)s" % dict(num=self.num, ext=self.config["Plots"]["format"],**self.config["Dirs"]))
        plt.clf()

        
    
    def plot_trace(self):
        """Two plots showing the trace data for this lenslet.
        
        1. Delta Lambda per Pixel vs. Wavelenth ``%(Partials)s/Instrument-%(num)04d-DeltaWL%(ext)s``
        
        2. Sampling Resolution vs. Wavelength ``%(Partials)s/Instrument-%(num)04d-Resolution%(ext)s``
        
        **Variables which are used**:
        
        :var tfl: flux of each pixel in counts
        :var twl: wavelength of each pixel in meters
        :var tdw: delta wavelength covered by each pixel in meters
        :var trs: sampling resolution of each pixel
        
        """
        assert self.traced
        plt.clf()
        plt.plot(self.twl*1e6,self.tdw*1e6,"g.")
        plt.title("$\Delta\lambda$ for each pixel (%d)" % self.num)
        plt.xlabel("Wavelength ($\mu m$)")
        plt.ylabel("$\Delta\lambda$ per pixel")
        plt.savefig("%(Partials)s/Instrument-%(num)04d-DeltaWL%(ext)s" % dict(num=self.num, ext=self.config["Plots"]["format"],**self.config["Dirs"]))
        plt.clf()
        self.log.debug(npArrayInfo(self.trs,"Trace RS"))
        plt.semilogy(self.twl*1e6,self.trs,"g.")
        plt.title("$R = \\frac{\lambda}{\Delta\lambda}$ for each pixel (%d)" % self.num)
        plt.xlabel("Wavelength ($\mu m$)")
        plt.ylabel("Resolution $R = \\frac{\lambda}{\Delta\lambda}$ per pixel")
        plt.savefig("%(Partials)s/Instrument-%(num)04d-Resolution%(ext)s" % dict(num=self.num, ext=self.config["Plots"]["format"],**self.config["Dirs"]))
        plt.clf()
    
    def bin(self,array,factor):
        """Bins an array by the given factor in both directions.
        
        :param array: array to be binned
        :param factor: binning factor
        :returns: array binned
        
        """
    
        finalShape = tuple((np.array(array.shape) / factor).astype(np.int))
        Aout = np.zeros(finalShape)
    
        for i in range(factor):
            Ai = array[i::factor,i::factor]
            Aout += Ai
    
        return Aout
        
    
    def rotate(self,point,angle,origin=None):
        """Rotate a given point by the provided angle around the origin given.
        
        :param point: The point in [x,y]
        :param angle: The angle in radians
        :param origin: The rotation origin for the point
        :returns: Rotated point [x,y]
        
        """
        if origin == None:
            origin = np.array([0,0])
        pA = np.matrix((point - origin))
        R = np.matrix([[np.cos(angle),-1*np.sin(angle)],[np.sin(angle),np.cos(angle)]]).T
        pB = pA * R
        return np.array(pB + origin)[0]
    
    def make_hexagon(self):
        """Generates a hexagonal polygon for this lenslet, to be used for area calculations.
        
        Generates a single point at the appropriate radius from the center. This point is then rotated by the hexagon rotation value. Next, five more points are generated by rotating the original point successive times by pi/3.
        
        :var shape: shapely polygon representing this hexagon.
        """
        radius = self.config["Instrument"]["Lenslets"]["radius"]
        angle = np.pi/3.0
        rotation =  self.config["Instrument"]["Lenslets"]["rotation"] * (np.pi/180.0) #Absolute angle of rotation for hexagons
        A = self.rotate(self.ps + np.array([radius,0]),rotation,self.ps)
        points = [A]
        for i in range(5):
            A = self.rotate(A,angle,self.ps)
            points.append(A)
        self.shape = sh.geometry.Polygon(tuple(points))
        
    def show_geometry(self,color='#cccc00',label=False):
        """Show the source geometry on the current plot. Geometry is shown filled with the color provided, and bordered by a grey box."""
        x, y = self.shape.exterior.xy
        plt.fill(x, y, color=color, aa=True) 
        plt.plot(x, y, color='#666666', aa=True, lw=0.25)
    
    def setup_crosstalk(self,n):
        """Set up the crosstalk matrix for this lenslet and the flat base spectrum for this lenselt. This setup depends on the intended size of the crosstalk matrix.
        
        :param n: Size of the crosstalk matrix (number of source pixels)"""
        self.pixelValues = np.zeros((n))
    
    def find_crosstalk(self,pixel):
        """Find the crosstalk overlap with a given pixel. The crosstalk is defined as the fraction of that pixels light which will end up in this lenslet.
        
        :param pixel: The pixel for calculation. :class:`SourcePixel`
        
        If the two shapes are disjoint, nothing is done. If they are not, we calculate the percentage of the pixel which ends up in this hexagon. This overlap factor is then saved in the crosstalk matrix for both the pixel and the lenslet. Finally, the pixel is scaled by the overlap and added to this lenslet's spectrum."""
        if self.shape.disjoint(pixel.shape):
            return
        overlap = (self.shape.intersection(pixel.shape).area) / pixel.shape.area
        self.spectrum += pixel * overlap
        self.pixelValues[pixel.idx] = overlap
        pixel.pixelValues[self.idx] = overlap
        
    def find_normalized_overlap(self):
        """Saves a normalized overlap matrix. The normalized overlap matrix is required for requesting colors from the color-map for matrix plotting."""
        matrix = self.pixelValues
        self.Norm_overlaps = Normalize(matrix)

        
class SourcePixel(InterpolatedSpectrum):
    """Source Pixels are objects which handle the source shape and size for resampling"""
    def __init__(self,x,y,config,num,**kwargs):
        super(SourcePixel, self).__init__(**kwargs)
        self.x = x
        self.y = y
        self.config = config
        self.ps = np.array([x,y])
        self.num = num
        self.method = self.resolve_and_integrate
    
    def make_pixel_square(self):
        """Make a specific pixel square"""
        radius = self.config["Source"]["PXSize"]["mm"]
        rotation = self.config["Source"]["Rotation"]
        angle = np.pi/2.0
        A = self.rotate(self.ps+np.array([0,radius]),rotation,self.ps)
        points = [A]
        for i in range(3):
            A = self.rotate(A,angle,self.ps)
            points.append(A)
        self.shape = sh.geometry.Polygon(tuple(points))
        
    def rotate(self,point,angle,origin=None):
        """Rotate a given point by the provided angle around the origin given"""
        if origin == None:
            origin = np.array([0,0])
        pA = np.matrix((point - origin))
        R = np.matrix([[np.cos(angle),-1*np.sin(angle)],[np.sin(angle),np.cos(angle)]]).T
        pB = pA * R
        return np.array(pB + origin)[0]
        
    def show_geometry(self,color='#cccc00'):
        """Show the source geometry"""
        x, y = self.shape.exterior.xy
        plt.fill(x, y, color=color, aa=True) 
        plt.plot(x, y, color='#666600', aa=True, lw=0.25)
        
    def setup_crosstalk(self,n):
        """docstring for setup_crosstalk"""
        self.pixelValues = np.zeros((n))
            
    def get_color(self,idx):
        """docstring for get_color"""
        self.vals = Normalize()(self.pixelValues)
        return cm.jet(self.vals[idx])
    
