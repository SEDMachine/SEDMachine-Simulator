
import numpy as np
from scipy.interpolate import interp1d

# electron volt, reflectivity
reflectivity_ev = np.array([
	[0.2, 0.9873],
	[0.25, .9858],
	[0.3, .9844],
	[0.4, .9826],
	[0.5, .9817],
	[0.6, .9806],
	[0.7, .9794],
	[0.8, .9778],
	[0.9, .9749],
	[1.0, .9697],
	[1.1, .9630],
	[1.2, .9521],
	[1.3, .9318],
	[1.4, .8852],
	[1.5, .8678],
	[1.6, .8794],
	[1.7, .8972],
	[1.8, .9069],
	[1.9, .9116],
	[2.0, .9148],
	[2.2, .9200],
	[2.4, .9228],
	[2.6, .9238],
	[2.8, .9242],
	[3.0, .9241],
	[3.2, .9243],
	[3.4, .9245],
	[3.6, .9246],
	[3.8, .9247],
	[4.0, .9248],
	[4.2, .9248],
	[4.4, .9249],
	[4.6, .9249],
	[4.8, .9249],
	[5.0, .9244],
	[6.0, .9257] ])

RD = reflectivity_ev
hc = 12398. #ev Ang
RD[:,0] = hc/RD[:,0]

s = np.argsort(RD[:,0])
interpfun = interp1d(RD[s,0], RD[s,1])

def reflectivity(lam):
	'''Reflectivity of Al as function of wavelength (Ang)'''

	return interpfun(lam)/interpfun(6700)*0.86



