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

test = SEDImage()
model = SEDSystem()
density = 5
size = model.npix

model.loadOpticsData("Data/lenslet_xy_13-10-2011.txt","Data/dispersion_12-10-2011.txt")
test.generate_blank((size,size))
points, intpoints, wl, deltawl, density = model.get_wavelengths(2128,density)
radiance = bbody(wl*1e-6) * 1e-6
flux = radiance[1,:-1] * deltawl

x,y = intpoints[:-1].T.astype(np.int)

xint,yint = points.T.astype(np.int)

padding = 10

x -= np.min(x)
y -= np.min(y)

xdist = np.max(x)-np.min(x)
ydist = np.max(y)-np.min(y)

xdist += (5 - xdist % 5)
ydist += (5 - ydist % 5)

x += padding * density
y += padding * density

img = np.zeros((xdist+2*padding*density,ydist+2*padding*density))

print img.shape

img[x,y] = flux

img2 = model.blur_image(img,12)

print img2.shape

small = bin(img2,density).astype(np.int16)


corner = np.array([ xint[np.argmax(x)], yint[np.argmin(y)]]) - np.array([padding,padding])

xstart = corner[0]
xend = xstart + small.shape[0]
ystart = corner[1]
yend = ystart+small.shape[1]

data = test.data()

print xstart,xend,ystart,yend,data.shape

data[xstart:xend,ystart:yend] += small

plt.figure()
plt.imshow(small)
plt.figure()
plt.imshow(img2)
plt.figure()
plt.imshow(img)
plt.show()

# test.place(2130,bbody,"Test",model,1)

# test.write("Image2.fits")


# test.place(2129,bbody,"HMMM",model)
# 
# print test.used
# 
# test.remove("Blank")
# test.crop(model.center[0],model.center[1],1024,1024)
# test.remove("HMMM")
# test.save(blur_image(test.data(),3),"Blur")
# os.remove("RotTest.fits")
# test.write("RotTest.fits")


