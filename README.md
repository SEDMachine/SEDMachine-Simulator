This is an introduction. A [full documentation][Docs] is on Github Pages

# Introduction to the SEDMachine Simulator

The [SEDMachine Simulator][Docs] is designed to simulate many aspects of the design and operation of the SEDMachine instrument. As such, it endeavors to be an accurate simulation of the optical processes in the SEDMachine and to faithfully reproduce data which might be created by the SEDMachine. This program can simulate the different types of data images that the SEDMachine is likely to produce. Images can be simulated which are based on data-cubes, simple sources (with a fixed geometry inserted), or uniform sources. Potential uniform sources include sky spectra, flat lamps (similar to what a quartz lamp would produce), calibration lamps (based on a line list) and arbitrary uniform sources (specified in a data file, with R>1000). The simulator works by generating large FITS images, and doing spectral resolution and flux calibration math. It can take a little while to run each simulation (there are more than 5000 active lenslets in the SED machine!), so be patient, and test your simulation on a small number of lenslets first to ensure that it works.

## Organization

The program is contained within the `SEDMachine` module. This module has two primary components: the Simulator, and the supporting classes. The `SEDMachine.Simulator` module has a single class, `SEDMachine.Simulator.SEDSimulator` which implements the control functions for the program. When you call `SEDMsim`, you are using this class's `SEDMachine.Simulator.SEDSimulator.run()` method. The program relies on the [AstroObject][] module, which I've (Alex Rudy) developed to make lots of astronomy tasks a lot more "object-focused". [AstroObject][] can be found online at [Github][].

There are also some support tools, including a basic data setup for the latest SEDMachine simulations. These support tools are kept in `SEDTools`.

This program is stored online at [Github][Repo], and the [documentation][Docs] is available there as well.

## What this program simulates

This program uses a ray-trace output from an optical model of the SEDMachine to simulate data as it will appear when read from the SEDMachine's CCDs. The simulator is only focused on simulating the work of the integral-field-unit and does not simulate or produce rainbow-camera data. The simulator accounts for the resolution of the instrument, the spectral profile, and the point-spread function of the system.

The simulation makes a few assumptions and uses a few "fudge" factors. These assumptions are documented in the `Model` which gives an overview of the modelling system, step-by-step.

## What this program does not do

This program does not handle ray-tracing, in any form. The ray tracing (and as such, the layout and physics of the system) are modeled separately and then interpreted by this simulator. The simulator also does not do rigorous scattered light calculations, nor does it manage the subtleties of image geometry.

## Very Basic Tutorial

This is the bare minimum required to just run a simulation. For a more detailed tutorial, see [Tutorial][].

Assuming you have [Matplotlib][], [Numpy][] and [Scipy][] installed, run:
	
	$ sudo python setup.py install
	$ cd working/directory/path/
	$ SEDMsetup *all *config-files
	$ SEDMsim *all
	
to perform a basic simulation. Simulations can take quite a long time to run. As such, you might want to try:
	
	$ SEDMsim -T *all
	
to perform the same simulation, but using only 50 lenslets.

## Building Documentation

If you can't find the [documentation][Docs] online (or it has moved), then:

1. Ensure you have Sphinx installed.
2. Navigate to the `Docs/` directory in the source tree you downloaded
3. Run `make html`
4. Open `Docs/build/html/index.html`

You should then have a complete usable HTML copy of the documentation.

[Tutorial]: http://alexrudy.github.com/SEDMachine-Simulator/
[Repo]: http://github.com/alexrudy/SEDMachine-Simulator/
[Docs]: http://alexrudy.github.com/SEDMachine-Simulator/
[Matplotlib]: http://matplotlib.sourceforge.net/
[Numpy]: http://numpy.scipy.org/
[Scipy]: http://scipy.org/
[MacPorts]: http://macports.org/
[AstroObject]: http://github.com/alexrudy/AstroObject/
[GitHub]: http://github.com/alexrudy/AstroObject/
