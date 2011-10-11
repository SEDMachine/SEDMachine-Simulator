#!/usr/bin/env python
# 
#  plotWL.py
#  SEDMachine
#  Version: 0.0
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 


import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pyfits
from scipy.interpolate import splprep, splev, splrep, interp1d
print("[plotWL] IMPORT COMPLETE")


showResiduals = True
showImages = True
showInterpolation = True
showSpectra = True
figureNumber = 1
# Read raw data
wavelengthTable = np.genfromtxt("Rtest.dat")
x,y = wavelengthTable[:,0],wavelengthTable[:,1]
xdense = np.arange(min(x),max(x),(max(x)-min(x))/1000000)
ydense = np.arange(min(y),max(y),(max(y)-min(y))/1000000)


# Spline interpolation
k = 3 # Order
s = 3.0 # Smoothness Parameter

# Compute the appropriate spline function
(t,c,k) = splrep(y,x)

# Evaluate the spline function at each data point
xint = splev(ydense,(t,c,k))

fspline = lambda ay : splev(ay,(t,c,k))


# Compute a linear interpolation function
flinear = interp1d(y,x)

#FITS File Generation
ImageSize = (500,20)
# An empty 500x20 Image array
SplineInterpolatedImage = np.zeros(ImageSize)
LinearInterpolatedImage = np.zeros(ImageSize)
NoInterpolationImage = np.zeros(ImageSize)

# The spectra are 1-D arrays of data which might be plotted later, instead of images.
# They are limited to the lenght of data, and have no padding around the edges.
SplineInterpolatedSpectra = []
LinearInterpolatedSpectra = []
NoInterpolationSpectra = []

# Iterate across all of the points which we are interested in.
for j in range(9,13):
    for i in range(ImageSize[0]):
        # How far (in Microns) have we moved at each pixel position
        distance = (i-120) * 13.5 / 1000
        
        # If this movement is outside the understood conversion function, ignore it.
        if distance > y[-1] or distance < y[0]:
            NoInterpolationImage[i,j] = 0
            LinearInterpolatedImage[i,j] = 0
            SplineInterpolatedImage[i,j] = 0
        # Else, we should populate the data
        else:
            # Find teh nearest position for this distance
            # This is a Floor interpolation function
            
            idx=(np.abs(y-distance)).argmin()
            NoInterpolationImage[i,j] = x[idx]
            
            # This uses the previously generated linear interpolation
            LinearInterpolatedImage[i,j] = flinear(distance)
            
            # This uses a generated spline interpolation
            SplineInterpolatedImage[i,j] = fspline(distance)
            
for i in range(ImageSize[0]):
    distance = (i-120) * 13.5 / 1000
    if distance > y[-1] or distance < y[0]:
        pass
    else:
        idx=(np.abs(y-distance)).argmin()
        NoInterpolationSpectra.append(x[idx])
        LinearInterpolatedSpectra.append(flinear(distance))
        SplineInterpolatedSpectra.append(fspline(distance))

# Display the prospective FITS image spectra
if showImages:
    plt.figure(figureNumber)
    figureNumber += 1
    plt.subplot(3,1,1)
    plt.title("Not Interpolated Spectra: Floor")
    plt.imshow(NoInterpolationImage.T)
    plt.subplot(3,1,2)
    plt.title("Interpolated Spectra: Linear")
    plt.imshow(LinearInterpolatedImage.T)
    plt.subplot(3,1,3)
    plt.title("Interpolated Spectra: Spline")
    plt.imshow(SplineInterpolatedImage.T)

# Display 1-D spectra
if showSpectra:
    plt.figure(figureNumber)
    figureNumber +=1
    d = (np.arange(len(NoInterpolationSpectra)) - 120) * 13.5 / 1000
    plt.title("Interpolated Spectra")
    plt.plot(d,NoInterpolationSpectra,'k-',label="No Interpolation")
    plt.plot(d,LinearInterpolatedSpectra,'g-')
    plt.plot(d,SplineInterpolatedSpectra,'b-')
    plt.xlabel("Wavelength")
    plt.ylabel("Distance along CCD")
    
# Display the interpolation residuals
if showResiduals:
    plt.figure(figureNumber)
    figureNumber +=1
    if showInterpolation:
        plt.subplot(3,1,2)
    else:
        plt.subplot(2,1,2)
    y = np.array(SplineInterpolatedSpectra)-np.array(NoInterpolationSpectra)
    y2 = np.array(LinearInterpolatedSpectra)-np.array(NoInterpolationSpectra)
    d = (np.arange(len(NoInterpolationSpectra)) - 120) * 13.5 / 1000
    plt.plot(d,y,'r.')
    plt.plot(d,y2,'g.')
    ys = y + y2
    plt.axis([min(d)-(max(d)-min(d))*0.05, max(d)*1.05, min(ys)-(max(ys)-min(ys))*0.05, max(ys)*1.05])
    # plt.xlabel("Wavelength")
    plt.ylabel("Residulas")
    plt.title("Difference between interpolation and Floor in microns of wavelength")
    if showInterpolation:
        plt.subplot(3,1,3)
    else:
        plt.subplot(2,1,1)
    y = np.array(SplineInterpolatedSpectra)-np.array(LinearInterpolatedSpectra)
    plt.plot(d,y,'r.')
    plt.axis([min(d)-(max(d)-min(d))*0.05, max(d)*1.05, min(y)-(max(y)-min(y))*0.05, max(y)*1.05])
    plt.xlabel("Distance along CCD")
    plt.ylabel("Residulas")
    plt.title("Difference between Spline and Linear in microns of wavelength")
        
    if showInterpolation:
        plt.subplot(3,1,1)
    
# Display the interpolation
if showInterpolation:
    if not showResiduals:
        plt.figure(figureNumber)
        figureNumber +=1
    x,y = wavelengthTable[:,0],wavelengthTable[:,1]
    plt.plot(x,y,'k.')
    plt.plot(xint,ydense,'b-')
    plt.plot(flinear(ydense),ydense,'g-')
    plt.axis([min(x)-(max(x)-min(x))*0.05, max(x)*1.05, min(y)-(max(y)-min(y))*0.05, max(y)*1.05])
    plt.xlabel("Wavelength")
    plt.ylabel("Distance along CCD")
    plt.title("Wavelength dispersion function on SED Machine CCD")

# Difference between Spline and Data at known points
xtest = splev(y,(t,c,k))

plt.figure(5)
plt.plot(y,xtest-x)
plt.xlabel("Distance along CCD")
plt.ylabel("Difference in Wavelength Spline-Given")
plt.title("Spline Residual from Calibration")

plt.show()
