import numpy as np
import scipy as sp
import pylab as pl
import scipy.interpolate, scipy.integrate

from matplotlib.backends.backend_pdf import PdfPages

import datetime

pl.rcParams.update({"legend.fontsize": 8})

def fraction_of_object(diam_as, seeing_FWHM_as):
	rad_as = diam_as/2.0
	fun = scipy.special.hyp2f1
	gamma = scipy.special.gamma
	FWHM = seeing_FWHM_as
	beta = 2.3

	coeff = 4*np.sqrt(2**(1/beta)-1) * gamma(beta)/(FWHM*np.sqrt(np.pi)*gamma(beta-0.5))
	p4 = -4* (2**(1/beta)-1)*rad_as**2/FWHM**2

	return coeff * rad_as * fun(0.5, beta, 1.5, p4)

def abmag_to_flambda(AB # Ab magnitude
	, lam): # Wavelength in angstrom

	c = 2.9979e18 # angstrom/s

	# return erg/s/cm^2/ang
	return 10**(-(AB + 48.6)/2.5)*c/lam**2

def vegamag_to_flambda(m, band):
	if band == 'r':
		pass

pl.ion()

hc = 1.98644521e-8 # erg angstrom


# See derivation on pg 83 of SED NB 1 (20 July 2011)
moon_phase = np.array([0., 0.08, 0.16, 0.24, 0.32, 0.40, 0.50])
moon_g = np.array([2e-17, 2.1e-17, 2.15e-17, 2.3e-17, 5.3e-17, 1.7e-16, 3.2e-16])
moon_r = np.array([2.3e-17,2.3e-17,2.3e-17,3.3e-17,3.5e-17,8.3e-17,1.3e-16])
moon_i = np.array([2.8e-17,3.0e-17,3.0e-17,3.3e-17,3.8e-17,7.0e-17,9.0e-17])

sky_ls = (4868., 6290., 7706., 10000)


moon_funs = []
for i in xrange(len(moon_phase)):
	gm = moon_g[i]-moon_g[0]
	rm = moon_r[i]-moon_r[0]
	im = moon_i[i]-moon_i[0]
	zm = im

	ff= np.poly1d(np.polyfit(sky_ls, np.array([gm, rm, im, zm]), 2))

	moon_funs.append(ff)

# From Turnrose PASP 86 (1974)
# Wavelength [A],  10^-18 erg / sec / cm^2 / Angstrom / Arcsecond^2
skyspec = [
	(3180, 4.30),
	(3260, 5.18),
	(3340, 6.13),
	(3420, 4.75),
	(3500, 4.86),
	(3580, 5.29),
	(3660, 7.24),
	(3740, 4.75),
	(3820, 4.43),
	(3900, 3.45),
	(3980, 4.31),
	(4060, 8.58),
	(4140, 6.09),
	(4220, 5.83),
	(4300, 5.39),
	(4380, 11.40),
	(4460, 6.25),
	(4540, 6.38),
	(4620, 6.16),
	(4700, 6.27),
	(4780, 6.14),
	(4860, 6.45),
	(4940, 6.24),
	(5020, 5.60),
	(5100, 5.80),
	(5180, 6.37),
	(5260, 6.26),
	(5340, 6.56),
	(5420, 7.85),
	(5500, 11.00),
	(5580, 25.40),
	(5660, 7.78),
	(5740, 9.70),
	(5760, 9.43),
	(5920, 11.40),
	(6080, 7.89),
	(6240, 13.00),
	(6400, 9.60),
	(6560, 8.36),
	(6720, 6.67),
	(6880, 9.73),
	(7040, 7.11),
	(7200, 9.53),
	(7360, 13.80),
	(7520, 10.70),
	(7680, 13.20),
	(7840, 23.60),
	(8000, 16.60),
	(8160, 5.54),
	(8320, 22.70),
	(8480, 19.30),
	(8640, 20.10),
	(8800, 36.10),
	(8960, 28.30),
	(9120, 8.22),
	(9280, 21.40),
	(9440, 32.40),
	(9600, 15.80),
	(9760, 26.30),
	(9920, 66.00),
	(10080, 68.30),
	(10240, 99.60),
	(10400, 87.10),
	(10560, 25.80),
	(10720, 64.30),
	(10880, 134.00)
]


pl.ion()

skyspec = np.array(skyspec)
skyspec[:,1] *= 1e-18 * 3 

def sky_function(PHASE):
	sky = skyspec[:,1] + moon_funs[PHASE](skyspec[:,0])
	sky /= hc/skyspec[:,0] # #/s/cm^2/ang/as^2
	skyf = sp.interpolate.interp1d(skyspec[:,0], sky)

	return skyf

import detectors
reload(detectors)
cameras = detectors.cameras

import atmosphere
reload(atmosphere)
ext = atmosphere.ext


# Source spectrum in units of 
# #/s/Angstrom/cm^2
# AB Mag at 6000
def sourcefun(ABMAG):

	return lambda x: abmag_to_flambda(ABMAG, x) / (hc/x)

thpts = np.load("thpt.npy")[0]
qeprism_pi = sp.interpolate.interp1d(thpts["lambda"], thpts["thpt-prism-PI"])
qeprism_andor = sp.interpolate.interp1d(thpts["lambda"], thpts["thpt-prism-Andor"])
qegrating = sp.interpolate.interp1d(thpts["lambda"], thpts["thpt-grating"])

def num_photons(fun, exptime, area, lam1, lam2):
	num_photon = sp.integrate.quad(fun, lam1, lam2, limit=150)[0] * exptime * area 

	return num_photon


def calculate(config, observation, figure):
	ls = []
	nsky_photon = []
	nsource_photon = []

	camera = cameras[config["camera_name"]]
	if config["camera_name"] == "PI":
		qe = qeprism_pi
	if config["camera_name"] == "Andor":
		qe = qeprism_andor
	if config["camera_name"] == "PIfast":
		qe = qeprism_pi

	if config["disperser"] == "grating":
		qe = qegrating

	RN = camera["RN"]
	DC = camera["DC"]
	R = config["R"]
	exptime = observation["exptime"] - camera["readtime"]*observation["n_exp"]
	skyf = sky_function(observation["PHASE_IX"])
	source = sourcefun(observation["source_AB"])
	telarea = observation["telarea_cm2"]

	startl = 8480
	endl = startl + startl/R

	spaxel_pitch_as = observation["spaxel_pitch_as"]
	seeing_fwhm = observation["seeing_fwhm_as"]
	seeing_area = (config["extract_diameter_spaxels"]*spaxel_pitch_as)**2 # as^2

	airmass = observation["airmass"]

	TOTAL_RN = np.sqrt(RN**2 * config["pix_per_res"] * 
		config["extract_pixs"] * config["extract_diameter_spaxels"]**2 * 
		observation["n_exp"] / config["GAIN"])
	TOTAL_DARK = DC * exptime

	# "slit loss" is 1-fraction
	fraction = fraction_of_object(spaxel_pitch_as*config["extract_diameter_spaxels"], 
		seeing_fwhm)

	while startl < 8680:
		l = (startl+endl)/2
		nsky_photon.append(num_photons(skyf, exptime, telarea, startl, endl) * \
			qe(l) * seeing_area)
		nsource_photon.append(num_photons(source, exptime, telarea, startl, endl) * \
			qe(l) * 10**(-ext(l)*airmass/2.5) * fraction)
		ls.append(l)

		startl = endl
		endl = startl + startl/R

	(ls, nsky_photon, nsource_photon) = map(np.array, (ls, nsky_photon, nsource_photon))
	shot_noise = np.sqrt(nsky_photon + nsource_photon + TOTAL_RN**2 + TOTAL_DARK)

	pl.figure(figure)
	pl.semilogy(ls, skyf(ls))
	pl.semilogy(ls, source(ls))

	feature_width = 20. # Angstrom
	instrument_width = 8500./R
	n_res = feature_width/R

	feature_photons = np.median(nsource_photon) * n_res
	feature_shot_noise = np.median(shot_noise) * np.sqrt(n_res)


	return (ls, nsky_photon, nsource_photon, shot_noise, TOTAL_DARK, 
		TOTAL_RN, fraction, qe(ls), 10**(-ext(ls) * airmass/2.5), 
		feature_photons, feature_shot_noise)

	
######################################
######################################
######################################
######################################

NPLOT = 10
for i in xrange(NPLOT+1): 
	pl.figure(i)
	pl.clf()

pl.figure(2)
pl.clf()
pl.ylim([1e1,1e6])
pl.grid(True, which='both')
pl.xlabel(u"$\lambda$ [$\AA$]")
pl.ylabel(u"$\gamma [e^-]$")
pl.title(u"Detected Counts from Sky (solid), Source (thick), \nRN$/G$ (straight), and Total Noise (dashed)")






legends= []
legends_two = []
for phiix in [0,1,2,3,4,5,6]:
	print "Phi: ", phiix

	config = {"GAIN": 2, "pix_per_res": 2.1, "extract_pixs": 4, "extract_diameter_spaxels": 1, "camera_name": "PI", 'disperser': 'prism', "R": 2000}
	observation = {"telarea_cm2": 45000, "exptime": 18000, "n_exp": 5, "airmass": 1.2, "PHASE_IX": phiix, "seeing_fwhm_as": 1.2, "spaxel_pitch_as": 0.65, "source_AB": 20}
	ls, nsky_photon, nsource_photon, shot_noise, dark, RN, fraction, qe, atm, feature_photons, feature_shot_noise = calculate(config, observation, 1)


	pl.figure(2)
	pl.semilogy(ls, nsky_photon)
	pl.semilogy(ls, shot_noise, '--')

	print "Total RN: %3.1f per res element" % RN
	print "Fraction of light captured %1.2f" % fraction

	pl.figure(3)
	pl.plot(ls, nsource_photon/shot_noise, linewidth=1)

	legends.append(u"$\phi$: %s" % moon_phase[observation["PHASE_IX"]])

pl.figure(1)
pl.grid(True, which='both')
pl.xlabel(u"$\lambda$ [$\AA$]")
pl.ylabel(u"$\gamma\ s^{-1} cm^{-2} \AA^{-1}$ (per $arcsec^{-2}$ for sky)")
pl.title("Photon Brightness from Sky and Object %2.1f" % observation["source_AB"])

pl.figure(2)
pl.semilogy(ls, nsource_photon, linewidth=3)
pl.semilogy(ls, np.ones(len(ls))*RN)

pl.figure(3)
pl.grid(True)
pl.ylim([0,15])
pl.xlabel(u"$\lambda$ [$\AA$]")
pl.ylabel("S/N per Resolution Element for AB %2.1f source" % observation["source_AB"])
pl.legend(legends, loc='upper left')
pl.title(u"%s: $m_{obj}=%2.1f$, D$_{extract}$=%3.1f, $f_{obs}$=%3.2f, $\Sigma RN$=%3.0f, $G$=%1.0f\nArea=%1.1e $cm^2$, exptime: %4i s, FWHM: %1.1f" % (config["camera_name"], observation["source_AB"], config["extract_diameter_spaxels"] * observation["spaxel_pitch_as"], fraction, RN, config["GAIN"], observation["telarea_cm2"], observation["exptime"], observation["seeing_fwhm_as"]))


pl.figure(4)
pl.clf()
pl.grid(True, which='both')
pl.xlabel(u"$\lambda$ [$\AA$]")
pl.ylabel("[fraction]")
pl.title("System Throughput: FWHM %2.1f\" / Extraction: %2.1f\" f: %1.2f" % (observation["seeing_fwhm_as"], config["extract_diameter_spaxels"] * observation["spaxel_pitch_as"], fraction))
pl.semilogy(ls, atm)
pl.plot(ls, qe)
pl.plot(ls, qe*atm)
pl.plot(ls, qe*atm*fraction)
pl.legend(["Atm", "Sys. QE", "Sys QE x Atm", "Sys QE x Atm x f"])
pl.ylim([.01, 1])


extract_widths = [1,2,3,4,5,6,7,8,9,10]


pl.figure(13)
pl.clf()
pl.figure(5)
pl.clf()
pl.grid(True)
pl.ylim([0,15])
pl.xlabel(u"$\lambda$ [$\AA$]")
pl.ylabel("S/N per Resolution Element")
legend_txt = [] 
subtitle_txt = u"$\Sigma$RNs: "
for extract in extract_widths:
	print "extract width: ",extract 

	observation["PHASE_IX"] = 1
	config["extract_diameter_spaxels"] = extract
	ls, nsky_photon, nsource_photon, shot_noise, dark, RN, fraction, qe, atm, feature_photons, feature_shot_noise = calculate(config, observation, 13)


	pl.figure(5)
	pl.plot(ls, nsource_photon/shot_noise)

	print "Total RN: %3.1f per res element" % RN
	print "Fraction of light captured %1.2f" % fraction
	legend_txt.append("%1.1f\": %0.2f" % (config["extract_diameter_spaxels"] * observation["spaxel_pitch_as"], fraction))

	subtitle_txt += "%3.0f " % RN

config["extract_diameter_spaxels"] = 3

now = datetime.datetime.now()
fname = now.strftime("VelMachine-%Y-%m-%d__%H-%M-%S") + ".pdf"
pl.figure(0)
pl.clf()
pos = .9
pl.title("Typical values: %s" % (__file__))
for key, value in config.items():
	pl.text(.2, pos, "%s: %s" % (key, value))
	pos -= .05

for key,value in observation.items():
	pl.text(.2, pos, "%s: %s" % (key, value))
	pos -= 0.05

camera = cameras[config["camera_name"]]
pl.text(.2, pos, "%s: %s" % ("DC", camera["DC"]))
pos -= 0.05
pl.text(.2, pos, "%s: %s" % ("RN", camera["RN"]))
pos -= 0.05

pos -= 0.05
pl.text(.2, pos, fname)

pl.figure(5)
pl.title(u"By Extraction Diameter -- %s: $m_{obj}=%2.1f$, $G$=%1.0f $\phi$=%0.2f\nArea=%1.1e $cm^2$, exptime: %4i s, FWHM: %1.1f" % (config["camera_name"], observation["source_AB"], config["GAIN"], moon_phase[observation["PHASE_IX"]], observation["telarea_cm2"], observation["exptime"], observation["seeing_fwhm_as"]))
pl.legend(legend_txt, loc='upper left')

pl.figure(6)
pl.xlabel(u"$R$")
pl.ylabel("Velocity Precision")
legend = []
Rs = np.array([ 10, 30, 100, 300, 1000, 3000])
fSNRs = []
widths = []
fCents = []

for R in Rs:
	config["R"] = R

	stellar_vel = 2.0

	width = np.sqrt(stellar_vel**2 + (3e5/R)**2)

	print R, width

	ls, nsky_photon, nsource_photon, shot_noise, dark, RN, fraction, qe, atm, feature_photons, feature_shot_noise = calculate(config, observation, 13)

	fSNR = feature_photons / feature_shot_noise
	fSNRs.append(fSNR)

	widths.append(width)


widths = np.array(widths)
fSNRs = np.array(fSNRs)

print widths

pl.figure(6)
pl.grid(True)
pl.semilogx(Rs, widths / fSNRs)



