#!/usr/bin/env python
# 
#  BlackBodyImage.py
#  Simulation Software
#  
#  Created by Alexander Rudy on 2011-10-21.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 

from SEDImage import *
import scipy.signal
import scipy.ndimage.interpolation
import os


bbody = BlackBodySpectrum(5000)
bbody1 = BlackBodySpectrum(100)
flat = FlatSpectrum(1e12)

image = SEDImage()
clip = ImageObject()
model = SEDSystem()

density = 5
size = model.npix
lenslet = 2129
spectrum = bbody
label = "BlackBody"
resolution_element = 2.4
psf_fwhm = 2.4
gain = 1e2

model.loadOpticsData("Data/lenslet_xy_13-10-2011.txt","Data/dispersion_12-10-2011.txt")
img,corner = model.get_dense_image(lenslet,spectrum,model,density)

clip.save(img,label + "-Dense")

img_circ = sp.signal.convolve(img,model.circular_kern(resolution_element*density/2.0),mode='same')

clip.save(img_circ,label+"-Dense-Circ")

img2 = model.blur_image(img_circ,psf_fwhm*density)

clip.save(img2,label+"-Dense-Blurry")

small = bin(img2/gain,density).astype(np.int16)

clip.save(small,label+"-Blurry")


plt.figure()
clip.show(label+"-Dense")
plt.savefig(label+"-Dense.png")
plt.figure()
clip.show(label+"-Blurry")
plt.savefig(label+"-Blurry.png")
plt.figure()
clip.show(label+"-Dense-Blurry")
plt.savefig(label+"-Dense-Blurry.png")
plt.figure()
clip.show(label + "-Dense-Circ")
plt.savefig(label+"-Dense-Circ.png")
plt.show()

clip.write("Lenslet%4d.fits" % lenslet,clobber=True)

