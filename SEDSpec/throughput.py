
import numpy as np
import pylab as pl
import aluminum as Al
import atmosphere
import detectors
import grating

from scipy.interpolate import interp1d

camera = {"Andor": { # Midband
		"DQEs": np.array([
		(2500, .05),
		(3500, .18),
		(4500, .75),
		(5000, .9),
		(5500, .92),
		(6500, .91),
		(7500, .79),
		(8500, .48),
		(9500, .13),
		(10500, .02),
		(11000, 0)
	]), 
	"RN": 2.9,
	"DC": 0.0004,
	"readtime":  82*3}}



l = np.arange(3200., 10501.)

# camcol: 14*2, Prism: 4, LA: 3 (scatter loss), Foreoptics: 6, Sp: 1, Detector: 2
nsurf = 14*2 + 4 + 3 + 6 + 1 + 2
nsurf_rc = 20

glass = l.copy()
glass_rc = l.copy()
glass[:] = .99**nsurf
glass_rc[:] = .99**nsurf_rc * .90


tel = Al.reflectivity(l)**2
airmass = 1.2
atm = 10**(-atmosphere.ext(l)*airmass/2.5)

qe = detectors.qe_function("Andor")
qe2 = detectors.qe_function("PI")


grat = grating.efficiency(l, 1, 5500, .9)/.99**3

fig = pl.figure(1)
pl.clf()
pl.ylim(.1,1)
pl.xlim([3700,9200])

pl.semilogy(l, tel)
pl.semilogy(l, atm)
pl.semilogy(l, tel*atm)
pl.semilogy(l, glass)
pl.semilogy(l, tel*atm*glass)
pl.semilogy(l, glass*qe2(l), lw=2)
pl.semilogy(l, tel*atm*glass*qe2(l), lw=4)
pl.xlabel(u"$\lambda$ [$\AA$]")
pl.ylabel("QE")
pl.yticks([.1,.2,.3,.4,.5,.6,.7,.8,.9,1], ("0.1", "0.2", "0.3", "0.4", "0.5",
    "0.6", "0.7", "0.8", "0.9", "1.0"))
pl.title("Predicted Throughput IFU")


pl.grid(True, which='both')

pl.legend(["Telescope", "Airmass=1.2", u"Tel$\\times$ Atm", "Optics",
    "On Detector", "IFU", u"Tel$\\times$ Atm$\\times$ IFU"], loc='lower center')


pl.savefig("predicted_throughput_8_nov_2011.pdf")

fig = pl.figure(2)
pl.clf()
pl.ylim(.01, 1)
pl.xlim([3200, 10000])

pl.semilogy(l, tel)
pl.semilogy(l, atm)
pl.semilogy(l, tel*atm)
pl.semilogy(l, glass_rc)
pl.semilogy(l, tel*atm*glass_rc)
pl.semilogy(l, glass_rc*qe2(l), lw=4)
pl.semilogy(l, tel*atm*glass_rc*qe2(l), lw=4)
pl.xlabel(u"$\lambda$ [$\AA$]")
pl.ylabel("QE")
pl.title("Predicted Throughput RC")


pl.grid(True, which='both')
pl.yticks([.1,.2,.3,.4,.5,.6,.7,.8,.9,1], ("0.1", "0.2", "0.3", "0.4", "0.5",
    "0.6", "0.7", "0.8", "0.9", "1.0"))

for WW in [3200, 3850, 4010, 5500, 5550, 6950, 6900, 8200]: pl.axvline(WW)

pl.legend(["Telescope", "Airmass=1.2", u"Tel$\\times$ Atm", "Optics", "On Detector", "RC", u"Tel$\\times$ Atm$\\times$ RC"], loc="lower center")
pl.savefig("predicted_throughput_RC_15_august_2011.pdf")

np.save("thpt", [{"lambda": l, "thpt-prism-PI": tel*atm*glass*qe2(l), "thpt-prism-Andor": tel*atm*glass*qe(l), "thpt-grating": tel*atm*glass*qe2(l)*grat, "thpt-rc-PI": tel*atm*glass_rc*qe2(l), "thpt-rc-Andor": tel*atm*glass_rc*qe(l)}])
