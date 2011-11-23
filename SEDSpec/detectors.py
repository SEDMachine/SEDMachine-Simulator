
import numpy as np
from scipy.interpolate import interp1d


cameras = {
	"PI": {"DQEs": np.array([
		(2000, 0),	
		(3000, 0.01),
		(3500, .20),
		(4000, .60),
		(4500, .82),
		(5000, .90),
		(5500, .93),
		(6000, .93),
		(7000, .93),
		(7500, .88),
		(8000, .73),
		(8500, .55),
		(9000, .33),
		(10000, .08),
		(10500, 0.02),
		(11000, 0)
	]),
	"RN" : 5.,
	"DC":  0.006,
	"readtime": 37},
	"Andor": { # Midband
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
	"RN": 4,
	"DC": 0.0004,
	"readtime":  82},
	"E2V": {"DQEs" : np.array([
		(3000, .1),
		(3500, .3),
		(4500, .8),
		(5500, .8),
		(6500, .78),
		(7500, .7),
		(8500, .4),
		(9500, .13),
		(10500, .02),
		(11000, 0)]),
	"RN": 3.3,
	"DC": 0.006,
	"readtime": 37},
}


cameras["PI-fast"] = cameras["PI"]
cameras["PI-fast"]["RN"] = 12
cameras["PI-fast"]["readtime"] = 2.265

cameras["Andor-fast"] = cameras["Andor"]
cameras["Andor-fast"]["RN"] = 11.7
cameras["Andor-fast"]["readtime"] = 1.398

def qe_function(camera_name):
	camera = cameras[camera_name]

	DQEs = camera["DQEs"]
	qe = interp1d(DQEs[:, 0], DQEs[:, 1], kind='cubic')

	return qe

