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
import os,time,sys

from multiprocessing import Pool, Value

bbody = BlackBodySpectrum(5000)
bbody1 = BlackBodySpectrum(100)

test = SEDImage()
model = SEDSystem()
model.loadOpticsData("Data/lenslet_xy_13-10-2011.txt","Data/dispersion_12-10-2011.txt")
test.generate_blank((model.npix,model.npix))
start = time.clock()
test.place_sed(2138,bbody,"BlackBody",model,5,3)
end = time.clock()
filename = "CCD_Bbody%05d.fits" % bbody.temperature
lenses = model.lenslets
# lenses = np.arange(50) + 2100

print "Placing a single spectrum took %3.4fs" % (end-start)
print "Placing all spectra should take about %6.4fs" % ((end-start) * len(model.lenslets))
print "That would be %d min" % ((end-start) * len(model.lenslets)/60)

last = 0
bar = arpytools.progressbar.ProgressBar()
print("Placing Spectra in %d Lenslets" % len(lenses))
bar.render(0,"%d" % 0)
total = float(len(lenses))
prog = Value('f',0.0)

def sed(i):
    """SED"""
    
    
    try:
        test.place_sed(i,bbody,"addOne%d" % i,model,5,3)
    except SEDLimits:
        LOG.warning("Couldn't place %d" % i)
    else:
        LOG.info("Placed %d" % i)
    finally:
        prog.value += 1.0
        bar.render(int(float(prog.value)/float(total) * 100.0),"%d" % prog.value)
        return

# start = time.clock()
# po = Pool()
# 
# 
# 
# po.map_async(sed,lenses)
# 
# po.close()
# po.join()

for i in lenses:
    sed(i)


end = time.clock()

stopw = end - start
totalt = (stopw / len(lenses)) * len(model.lenslets) / 60.0
nlenslets = len(model.lenslets)
ntestlens = len(lenses)
print "Placing %d spectra took %3.4fs,\n so all spectra should take about %3.4fmin (%d spectra)" % (ntestlens,stopw,totalt,nlenslets)
print np.max(test.data())
print "States Existing: ",test.list()
test.crop(model.center[0],model.center[1],1024,1024)
test.save(test.data().astype(np.int16).T,"Final")
test.keep("Final")
print "States to Save: ",test.list()

test.write("Test.fits",clobber=True)