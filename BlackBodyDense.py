#!/usr/bin/env python
# 
#  BlackBodyImage.py
#  Simulation Software
#  
#  Created by Alexander Rudy on 2011-10-21.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 

from SEDImage import *

bbody = BlackBodySpectrum(5000)
bbody1 = BlackBodySpectrum(100)

test = SEDImage()

test.loadOpticsData("Data/lenslet_xy_13-10-2011.txt","Data/dispersion_12-10-2011.txt")
test.over_dense_blank()

last = 0
disable_Console()
bar = arpytools.progressbar.ProgressBar()
print("Placing Spectra in All Lenslets")
bar.render(0,"%d" % 0)
total = float(max(test.lenslets))

for i in test.lenslets:
    bar.render(int(float(i)/total * 100),"%d" % i)
    try:
        test.over_dense_place(i,bbody,"addOne%d" % i)
        if last:
            test.remove("addOne%d" % last)
    except Exception:
        LOG.warning("Couldn't place %d" % i)
    else:
        LOG.info("Placed %d" % i)
        last = i
enable_Console()

test.over_dense_place(2129,bbody1,"twoBody")
test.remove("addOne%d" % last)
test.crop(test.center[0],test.center[1],1024,1024)
test.save(test.data().astype(np.int16),"Final")
test.remove("twoBody")
test.remove("Cropped")
test.remove("Blank-Dense")
print "States to Save: ",test.list()
try:
    os.remove("./TestBig.fits")
except:
    print("Couldn't Remove TestBig.fits")
test.write("TestBig.fits")