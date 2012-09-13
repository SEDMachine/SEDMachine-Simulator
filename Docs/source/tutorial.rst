.. _Tutorial:

Basic Tutorial
==============

This tutorial should walk you through the installation and basic operation of the :program:`SEDMsim`. For more details, see :ref:`SEDMsim`.

Installation
------------

This module uses the ``setup.py`` paradigm for installation and operation. The included setup.py file should attempt to resolve the dependency on `AstroObject`_. It may also try to resolve dependencies on `Matplotlib`_, `Numpy`_ and `Scipy`_. It is adivsable to install these three modules **before** you attempt to install this program. They are specific and picky modules which often require material to be compiled on your machine. A package manager like `APT`_ or `MacPorts`_ might be useful here.

1. Install ``Python 2.7``. This program runs only in version 2.7 for now.
2. Install `Matplotlib`_, `Numpy`_ and `Scipy`_ in your prefered package manager.
3. Run ``setup.py`` in the SEDM directory ::
	
	$ sudo python setup.py install
	

This should complete your installation. Common problems are caused by the inability to install `AstroObject`_ correctly. If this happens, download `AstroObject`_ from GitHub (use the master branch, or look at the latest tag for the most up-to-date version), change to the AstroObject base directory, and run ``setup.py`` there ::
	
	$ sudo python setup.py install
	
Then return to step 3 above. As much as I try to make the installation of `AstroObject`_ painless, sometimes it doesn't happen automatically with the `setup.py` file provided for SEDMachine.

Running :program:`SEDMsim`
--------------------------

The simulator has two basic components: The simulation, and the environment. Your environment should be a working directory, but **not** your source directory.

Setup of your environment is handled by the ``SEDMsetup`` command::
	
	$ SEDMsetup
	

This creates all of the necessary data and directory files for use with the program. Next, you will need a configuration file. A basic configuration file that works with the data provided can be generated using::
	
	$ SEDMsetup *config-file
	

If you change values in the configuration file (Such as the directory information), you may want to run ``SEDMsetup`` again to make sure you have the correct data files available.

To run the program, choose a stage to run (see :option:`{stages}` for a list of basic stages) and then run the SEDMsim command. The simplest command is ``*all`` which performs the default simulation ::
	
	$ SEDMsim *all
	

Run this command with the ``-T`` option to limit the number of lenslets evaluated ::
	
	$ SEDMsim -T *all
	

The operating directory of the simulation should contain the following sub-directories which handle input and output:

	- ``Caches/`` which will contain cache files. Cache files include the ``.fits`` files for each individual spectra, useful to examine spectra in over-sampled view.
	- ``Data/`` which will contain the required data files. Some data files are distributed with the system.
	- ``Images/`` which is the output directory. The final ``.fits`` image will be saved there.
	- ``Logs/`` which contains simulation logs.
	- ``Partials/`` which will contain items produced part-way through the simulator. Many files will be placed there, including plots.

Examining :program:`SEDMsim` output
-----------------------------------

The primary output will appear in the ``Images/`` directory. This directory will contain ``.fits`` files which are full ouput. Each simulation produces two ``.fits`` files. One file includes a ``-deep`` in the filename. Both will contain today's date in the filename. The ``-deep`` file is an extended fits file with frames for a variety of stages in the simulation process. The other file is just the final output. If you use ``ds9``, you can examine the extended frames in the fits file using ::
	
	$ ds9 file-deep.fits[{0,1,3}]
	
for frames 0, 1, and 3. The final image will always be in frame 0. Other frames will use the FITS keyword ``LABEL`` to describe their contents. These frames may also use the ``OBJECT`` keyword (if it in use) to describe their contents. You can also view all of the frames using the ``ds9`` ``Open Other...`` menu item.

If you made plots, or wish to look at results from the simulation which are not the final image, examine the ``Partials/`` directory. If you just ran the simulation normally, you can find the following files in the partials directory:

- ``Instrument-Offsets.dat``: Pixel offset values for subimages
- ``Lenslets-Raw.dat``: Information about each lenslet, in human-readable form
- ``center.dat``: Information about the central lenslet
- ``SEDMachine.config.yaml``: Configuration file used for the simulation.

If you ran the ``*plot`` stage, you will find a plethora of plots about the entire simulation. As plots are produced for each lenslet, it is advisable to limit the number of lenslets used for the ``*plot`` command, for example, run ::
	
	$ SEDMsim -S *plot
	
to only examine 5 lenslets.

Configuring :program:`SEDMsim`
------------------------------

The SEDMsim program is primarily configuraed by a YAML file. The full configuration is described in :ref:`SEDM_Configuration`. A simple configuration file is shown below::
	
	Instrument:
	  density: 5
	  eADU: 0.03802
	Observation:
	  airmass: 1
	  exposure: 1200
	  number: 3
	  Moon:
	    Phase: 6
	  Sky:
	    Atmosphere: Atmosph
	    Use: TurnroseSKY
	Source:
	  Filename: Data/SNIa.R1000.dat
	

Each directive is as follows:

- ``Instrument`` contains instrument specific directives
	- ``density`` is the level of oversampling for each spectrum. Higher numbers produce more accurate data, but are slower.
	- ``eADU`` is the conversion from photons to counts.
- ``Observation`` Contains data about the individual observation.
	- ``airmass`` for the target, controls extinction
	- ``exposure`` time measured in seconds
	- ``number`` is the number of individual exposures which contribute to the overall exposure
	- ``Moon`` contains information about the moon.
		- ``Phase`` is an index for the moon phase.
	- ``Sky`` Controls directives about the sky model used.
		- ``Atmosphere`` is the label for the atmosphere spectrum ``.fits`` file.
		- ``Use`` is the label for the sky spectrum ``.fits`` file.
- ``Source`` contains information about the source in use.
	- ``Filename`` is the filename for 1D source spectra


Placing this structure into a configuration file named ``SED.main.config.yaml`` will change the operational values for that spectrum.

Creating alternative data
-------------------------

There are few other data frames that could be useful with the :program:`SEDMsim` simulator. The simulator can create Flat, Bias, Dark, Calibration and Sky frames which all have universal source properties and are required for most data reduction pipelines.

Flat Frames
~~~~~~~~~~~

Flat frames are substituted using the ``*flat-source`` stage. This stage uses the configuration values ``Source: Value:`` to get a flat value. To build a flat source frame, run ::
	
	$ SEDMsim *flat-source *all
	

Bias and Dark Frames
~~~~~~~~~~~~~~~~~~~~

Simple bias and dark frames (using the observatin characteristics in the configuration) can be created using ::
	
	$ SEDMsim *add-noise *save
	

Calibration Lamp Frames
~~~~~~~~~~~~~~~~~~~~~~~

Calibration lamp frames are based on a line list and properties about that line list. The configuration is ::
	
    WLCal:
      List: Data/Lines.dat
      sigma: 1.0e-09
      value: 100000000.0
	
and they can be run with ::
	
	$ SEDMsim *line-source *all
	

Sky Line Frames
~~~~~~~~~~~~~~~

Frames with pure sky can be created using ::
	
	$ SEDMsim *sky-source *all
	


.. _Matplotlib: http://matplotlib.sourceforge.net/
.. _Numpy: http://numpy.scipy.org/
.. _Scipy: http://scipy.org/
.. _MacPorts: http://macports.org/
.. _APT: http://en.wikipedia.org/wiki/Advanced_Packaging_Tool
.. _AstroObject: http://github.com/alexrudy/AstroObject/
.. _GitHub: http://github.com/alexrudy/AstroObject/
