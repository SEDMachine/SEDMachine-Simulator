.. _Introduction:

Introduction to the SEDMachine Simulator
========================================

The SEDMachine Simulator is designed to simulate many aspects of the design and operation of the SEDMachine instrument. As such, it endeavors to be an accurate simulation of the optical processes in the SEDMachine and to faithfully reproduce data which might be created by the SEDMachine. This program can simulate the different types of data images that the SEDMachine is likely to produce. Images can be simulated which are based on data-cubes, simple sources (with a fixed geometry inserted), or uniform sources. Potential uniform sources include sky spectra, flat lamps (similar to what a quartz lamp would produce), calibration lamps (based on a line list) and arbitrary uniform sources (specified in a data file, with R>1000). The simulator works by generating large FITS images, and doing spectral resolution and flux calibration math. It can take a little while to run each simulation (there are more than 5000 active lenslets in the SED machine!), so be patient, and test your simulation on a small number of lenslets first to ensure that it works.

Organization
------------

The program is contained within the :mod:`SEDMachine` module. This module has two primary components: the Simulator, and the supporting classes. The :mod:`SEDMachine.Simulator` module has a single class, :class:`SEDMachine.Simulator.SEDSimulator` which implements the control functions for the program. When you call :program:`SEDMsim`, you are using this class's :meth:`SEDMachine.Simulator.SEDSimulator.run` method. The program relies on the :mod:`AstroObject` module, which I've (Alex Rudy) developed to make lots of astronomy tasks a lot more "object-focused". `AstroObject`_ can be found online at `Github`_

There are also some support tools, including a basic data setup for the latest SEDMachine simulations. These support tools are kept in :mod:`SEDTools`.

What this program simulates
---------------------------

This program uses a ray-trace output from an optical model of the SEDMachine to simulate data as it will appear when read from the SEDMachine's CCDs. The simulator is only focused on simulating the work of the integral-field-unit and does not simulate or produce rainbow-camera data. The simulator accounts for the resolution of the instrument, the spectral profile, and the point-spread function of the system.

The simulation makes a few assumptions and uses a few "fudge" factors. These assumptions are documented in :ref:`Model` which gives an overview of the modelling system, step-by-step.

What this program does not do
-----------------------------

This program does not handle ray-tracing, in any form. The ray tracing (and as such, the layout and physics of the system) are modeled separately and then interpreted by this simulator. The simulator also does not do rigorous scattered light calculations, nor does it manage the subtleties of image geometry.

Very Basic Tutorial
-------------------

This is the bare minimum required to just run a simulation. For a more detailed tutorial, see :ref:`Tutorial`.

Assuming you have `Matplotlib`_, `Numpy`_ and `Scipy`_ installed, run::
	
	$ sudo python setup.py install
	$ SEDMsetup *all *config-files
	$ SEDMsim *all
	
to perform a basic simulation. Simulations can take quite a long time to run. As such, you might want to try::
	
	$ SEDMsim -T *all
	
to perform the same simulation, but using only 50 lenslets.


.. _Matplotlib: http://matplotlib.sourceforge.net/
.. _Numpy: http://numpy.scipy.org/
.. _Scipy: http://scipy.org/
.. _MacPorts: http://macports.org/
.. _APT: http://en.wikipedia.org/wiki/Advanced_Packaging_Tool
.. _AstroObject: http://github.com/alexrudy/AstroObject/
.. _GitHub: http://github.com/alexrudy/AstroObject/