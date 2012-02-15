---
Title: Simulation Scripts for SED Machine
Author: Alexander Rudy

  Version 0.3.0

---

# Introduction
This Program contains SEDMachine Simulation and plotting scripts. They are written in Python. These scripts rely on the AstroObject modules that I wrote, which are bundled here. As well, they rely on matplotlib, numpy, scipy and pyfits. There is documentation for this program in the Docs folder.

# Release Notes

* 0.1.1
	- Initial Basic Release
* 0.1.2
	- Documentation of Command Line Interface
	- Update to SED.py script
	- Separation of Instrument functions into SEDInstrument.py
	- Creation of SEDSource module to handle source manipulation
	- Cleanup of Nick's spectrum simulations in SEDSpec
	- Basic Utlity additions
* 0.1.3
	- Use of setup.py for deployment and setup managment.
* 0.1.3p1,p2
	- Patch setup.py for AstroObject 0.2.8 rather than 0.2.6
	- Patch setup.py for Matplotlib 1.0 rather than 1.1
* 0.2.0a1
	- This alpha is simply a port to the new simulator framework. Only the main simulator has been ported at this point.
	- All modules will now use the framework-based logging system. This still produces slightly odd results.
	
