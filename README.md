---
Title: Simulation Scripts for SED Machine
Author: Alexander Rudy

  Version 0.3.2

---

# Introduction
This Program contains SEDMachine Simulation and plotting scripts. They are written in Python. These scripts rely on the AstroObject modules that I wrote, which are bundled here. As well, they rely on matplotlib, numpy, scipy and pyfits. There is documentation for this program in the Docs folder.

Documentation is avaliable on GitHub Pages: <http://alexrudy.github.com/SEDMachine-Simulator/>

# Installation

The program has a distribute-compatible setuputils script. To install the program, use:
	
	$ sudo python setup.py install
	
This will install two executable programs:

* `SEDMsetup` - A setup and data-file creation program
* `SEDMsim` - The simulator program

# Running the Program

Navigate your desired working directory, then run
	
	$ SEDMsetup
	
This will make working directories for your simulation.

To get an example configuration file which uses the included data, run
	
	$ SEDMsetup *config-file
		

Next run
	
	$ SEDMsim *all
	
to run the simulator.

For more information, see the documentation.

# Building Documentation

If you can't find the documentation online (or it has moved), then:

1. Ensure you have Sphinx installed.
2. Navigate to the `Docs/` directory in the source tree you downloaded
3. Run `make html`
4. Open `Docs/build/html/index.html`

You should then have a complete usable HTML copy of the documentation.