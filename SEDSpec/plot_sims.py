
import numpy as np
import scipy as sp
import pylab as pl

pl.figure(1)
pl.clf()
pl.grid(True)
pl.ylabel("Ratio of S/N")
pl.xlabel(u"$\lambda$ [$\AA$]")

pl.figure(2)
pl.clf()
pl.grid(True)
phases = ["0.00", "0.08", "0.16", "0.24", "0.32", "0.40", "0.50"]
pl.ylabel("Princeton PIXIS S/N per Resolution Element")
pl.title(u"Mag 21 source: 9 $arcsec^2$ area")
pl.xlabel(u"$\lambda$ [$\AA$]")
pl.ylim([0,16])

for phase in phases:

		dd = np.load("sn_%s_%s.npy" % ("Andor", phase))[0]
		sna = dd["source"] / dd["noise"]
		dd = np.load("sn_%s_%s.npy" % ("PI", phase))[0]
		snp = dd["source"] / dd["noise"]


	
		pl.figure(1)
		pl.plot(dd["lambda"], snp/sna)
		pl.ylim([.90,1/.90])

		pl.figure(2)
		pl.plot(dd["lambda"], snp)


pl.figure(2)
pl.legend(phases)
pl.savefig("s_n_9as.pdf")
pl.figure(1)
pl.legend(phases, loc="bottom right")
pl.title("Princeton/Andor SN")
pl.savefig("princeton_v_andor9_as.pdf")

