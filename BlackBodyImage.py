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
import os,time

from multiprocessing import Pool, Value


bbody = BlackBodySpectrum(5000)
bbody1 = BlackBodySpectrum(100)

test = SEDImage()
model = SEDSystem()
model.loadOpticsData("Data/lenslet_xy_13-10-2011.txt","Data/dispersion_12-10-2011.txt")
test.generate_blank((model.npix,model.npix))
start = time.clock()
test.place_sed(2138,bbody1,"",model,5,3)
end = time.clock()

print "Placing a single spectrum took %3.4fs" % (end-start)
print "Placing all spectra should take about %6.4fs" % ((end-start) * len(model.lenslets))
print "That would be %d min" % ((end-start) * len(model.lenslets)/60)

last = 0
bar = arpytools.progressbar.ProgressBar()
print("Placing Spectra in All Lenslets")
bar.render(0,"%d" % 0)
total = float(max(model.lenslets))
prog = Value('d',0.0)

def sed(i):
    """SED"""
    
    bar.render(int(prog.value/total * 100),"%d" % prog.value)
    try:
        test.place_sed(i,bbody,"addOne%d" % i,model,5,3)
    except SEDLimits:
        LOG.warning("Couldn't place %d" % i)
    else:
        LOG.info("Placed %d" % i)
    finally:
        prog.value += 1
        return

po = Pool()

po.map_async(sed,model.lenslets)

po.close()
po.join()
test.crop(model.center[0],model.center[1],1024,1024)
test.save(test.data().astype(np.int16).T,"Final")
test.remove("Cropped")
test.remove("Blank")
print "States to Save: ",test.list()
try:
    os.remove("Test.fits")
except:
    pass
test.write("Test.fits")