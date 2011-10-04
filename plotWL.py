#!/usr/bin/env python
# 
#  plotWL.py
#  SEDMachine
#  
#  Created by Alexander Rudy on 2011-10-03.
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
showInterpolation = False
showSpectra = True
figureNumber = 1
# Read raw data
wavelengthTable = np.genfromtxt("Rtest.dat")
x,y = wavelengthTable[:,0],wavelengthTable[:,1]
xdense = np.arange(min(x),max(x),(max(x)-min(x))/1000)
ydense = np.arange(min(y),max(y),(max(y)-min(y))/1000)


# Spline interpolation
k = 3 # Order
s = 3.0 # Smoothness Parameter

# Compute the appropriate spline function
(t,c,k) = splrep(y,x)

# Evaluate the spline function at each data point
xint = splev(ydense,(t,c,k))

# Compute a linear interpolation function
flinear = interp1d(y,x)

#FITS File Generation
ImageSize = (500,20)
# An empty 500x20 Image array
SplineInterpolatedImage = np.zeros(ImageSize)
LinearInterpolatedImage = np.zeros(ImageSize)
NoInterpolationImage = np.zeros(ImageSize)

SplineInterpolatedSpectra = []
LinearInterpolatedSpectra = []
NoInterpolationSpectra = []

for j in range(9,13):
    for i in range(ImageSize[0]):
        distance = (i-120) * 13.5 / 1000
        if distance > y[-1] or distance < y[0]:
            NoInterpolationImage[i,j] = 0
            LinearInterpolatedImage[i,j] = 0
            SplineInterpolatedImage[i,j] = 0
        else:
            idx=(np.abs(distance-y)).argmin()
            NoInterpolationImage[i,j] = x[idx]
            LinearInterpolatedImage[i,j] = flinear(distance)
            idx=(np.abs(distance-ydense)).argmin()
            SplineInterpolatedImage[i,j] = xint[idx]
            
for i in range(ImageSize[0]):
    distance = (i-120) * 13.5 / 1000
    if distance > y[-1] or distance < y[0]:
        pass
    else:
        idx=(np.abs(distance-y)).argmin()
        NoInterpolationSpectra.append(x[idx])
        LinearInterpolatedSpectra.append(flinear(distance))
        idx=(np.abs(distance-ydense)).argmin()
        SplineInterpolatedSpectra.append(xint[idx])

# Display the prospective FITS image spectra
if showImages:
    plt.figure(figureNumber)
    figureNumber += 1
    plt.subplot(3,1,1)
    plt.title("FITS Images of Interpolated Spectra")
    plt.imshow(NoInterpolationImage.T)
    plt.subplot(3,1,2)
    plt.imshow(LinearInterpolatedImage.T)
    plt.subplot(3,1,3)
    plt.imshow(SplineInterpolatedImage.T)

        
if showSpectra:
    plt.figure(figureNumber)
    figureNumber +=1
    d = (np.arange(len(NoInterpolationSpectra)) - 120) * 13.5 / 1000
    plt.plot(d,NoInterpolationSpectra,'k-')
    plt.plot(d,LinearInterpolatedSpectra,'g-')
    plt.plot(d,SplineInterpolatedSpectra,'b-')
    
# Display the interpolation residuals
if showResiduals:
    plt.figure(figureNumber)
    figureNumber +=1
    if showInterpolation:
        plt.subplot(2,1,2)
    y = np.array(LinearInterpolatedSpectra)-np.array(SplineInterpolatedSpectra)
    d = (np.arange(len(NoInterpolationSpectra)) - 120) * 13.5 / 1000
    plt.plot(d,y,'r-')
    plt.axis([min(d)-(max(d)-min(d))*0.05, max(x)*1.05, min(y)-(max(y)-min(y))*0.05, max(y)*1.05])
    plt.xlabel("Wavelength")
    plt.ylabel("Residulas")
    plt.title("Residual from Spline and Linear Fits")
    if showInterpolation:
        plt.subplot(2,1,1)
    
# Display the interpolation
if showInterpolation:
    if not showResiduals:
        plt.figure(figureNumber)
        figureNumber +=1
    plt.plot(x,y,'k.')
    plt.plot(xint,ydense,'b-')
    plt.plot(flinear(ydense),ydense,'g-')
    plt.axis([min(x)-(max(x)-min(x))*0.05, max(x)*1.05, min(y)-(max(y)-min(y))*0.05, max(y)*1.05])
    plt.xlabel("Wavelength")
    plt.ylabel("Distance along CCD")
    plt.title("Wavelength dispersion function on SED Machine CCD")

plt.show()
