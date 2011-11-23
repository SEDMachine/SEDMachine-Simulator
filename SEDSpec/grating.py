
import numpy as np

sin = np.sin
pi = np.pi

def efficiency(lam, order, blaze, peak):
	'''lam is wavelength in wavelength units, order is grating order, blaze is in wavelength units, peak is on blaze efficiency'''

	x = pi* (blaze/lam - order)

	return peak * (sin(x)/x)**2

