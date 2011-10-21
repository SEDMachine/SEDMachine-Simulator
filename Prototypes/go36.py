

# See notes on 12 Oct 2011 of SED Machine Notebook #2 pg 63
#


import numpy as np
import scipy as sp
import scipy.signal
import pylab as pl
import pyfits as pf

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = np.int(sizey)
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/np.float(size)+y**2/np.float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = sp.signal.convolve(im,g, mode='valid')
    return(improc)


try:
    f = open("../Data/lenslet_xy_13-10-2011.txt")
    lines = f.readlines()
    f.close()
except:
    print "Couldn't deal with file"

ix = []
p1 = []
p2 = []
xs = []
ys = []
lams = []
for line in lines[1:]:
    stp = map(float, line.rstrip().split())
    ix.append(stp[0])
    p1.append(stp[1])
    p2.append(stp[2])
    lams.append(stp[3])
    xs.append(stp[4])
    ys.append(stp[5])

[ix, p1, p2, lams, xs, ys] = map(np.array, [ix, p1, p2, lams, xs, ys])

print xs,ys
widthmm = 35.0

# It appears that xs,ys start as centered positions, (0-center of CCD), and are then adjusted to have a 0-0 starting point.
xs += widthmm/2
ys += widthmm/2

# Find the xs and ys that are not within 0.1 mm of the edge of the detector...
ok = (xs > 0.1) & (xs < widthmm-0.1) & (ys > 0.1) & (ys < widthmm-0.1)

# Slice all of the data arrays to account for this
ix = ix[ok]
p1 = p1[ok]
p2 = p2[ok]
lams = lams[ok]
xs = xs[ok]
ys = ys[ok]

# Calculate the number of pixels across the CCD
pixunit = 1 / (0.0135/1)
npix = np.round(widthmm * pixunit, 0)

# Spacing, in mm, between adjacent pixels, after rounding so that 35mm chip includes only whole pixels?
delt = widthmm / npix

# An empty image
img = np.zeros((npix, npix))

ls = np.arange(.37, .92, .001)

# Lenslet indicies
lenslets = np.unique(ix)
lengths = []

for lenslet in lenslets:
    # Indicies for selecting all data points which apply to a particular lenslet
    ok = ix == lenslet
    xout = np.round(xs[ok] * pixunit,0)
    yout = np.round(ys[ok] * pixunit,0)
    

    xo = yout.astype(np.int)
    yo = xout.astype(np.int)

    # Place the lenslet number into the positions for each lenslet
    img[xo[0],yo[0]-1] = lenslet
    img[xo[0],yo[0]+1] = lenslet
    img[xo[0]-1,yo[0]] = lenslet
    img[xo[0]+1,yo[0]] = lenslet

    
    
    if len(xo) < 3: 
        continue

    if np.any(xo == 0):
        continue

    if np.any(yo == 0):
        continue

    if np.any(np.abs(np.diff(yo)) > 30):
        continue

    lengths.append(xo[-1] - xo[0])

    fx = np.poly1d(np.polyfit(lams[ok], xo, 2))
    fy = np.poly1d(np.polyfit(lams[ok], yo, 2))
    if fx[1] > 3e3:
        print lenslet, fx[1]
        print xo
        print yo

    try: img[fx(ls).astype(np.int), fy(ls).astype(np.int)] = ls * 1000
    except: continue


cntix = np.argmin(p1**2 + p2**2)
center = (xs[cntix] * pixunit, ys[cntix] * pixunit)

img2 = img[center[0]-1024:center[0]+1024, center[1]-1024:center[1]+1024]
img2 = img2.astype(np.int16)


#bimg = blur_image(img, 3.0)
hdu = pf.PrimaryHDU(img2)

hdu.header.update("object", "Run 59g")
hdu.writeto("img40.fits",clobber=True)
