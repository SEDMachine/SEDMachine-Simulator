.. program:: SEDMsim
.. _SEDMsim:

:program:`SEDMsim` Program for running simulations
==================================================

The master simulator program is a command-line interface to the :meth:`AstroObject.AstroSimulator.Simulator.run` method. Below are the major command line components.

:program:`SEDMsim` Usage
------------------------

Usage Statement ::
	
	SEDMsim [ options ][ configuration ] {stages}
	
The program is actually agnostic to the order of arguments. Any argument may come in any position. As such, all arguments must be unique.


 .. option:: {stages}
	
	The stages option specifies individual stages for the program to run. You must specify at least one stage to run in the simulator. By default, two basic stages are provided, ``*all`` and ``*none``. The default simulation is performed by ``*all``. To test the simulator without running any stages (for example, to test :meth:`AstroObject.AstroSimulator.Simulator.registerFunction` functionality), us the ``*none`` stage to opertate without using any stages.
	
	Stages are called with either a ``*``, ``+`` or ``-`` character at the beginning. Their resepctive actions are shown below. All commands must include at least one macro. If you don't want any particular macro, use the ``*none`` macro.
	
	========= ============ ================================
	Character  Action      Description
	========= ============ ================================
	``*``     Include      To include a stage, use ``*stage``. This will also run the dependents for that stage.
	``-``     Exclude      To exclude a stage, use ``-stage``. This stage (and it's dependents) will be skipped.
	``+``     Include-only To include a stage, but not the dependents of that stage, use ``+stage``.
	========= ============ ================================
	
	.. Note ::
		Including an argument and excluding it simultaneously will have the effect of including it overall. So the following three commands are equivalent::
			
			$ SEDMsim +stage -stage *none
			$ SEDMsim -stage +stage *none
			$ SEDMsim +stage *none
			
		Also, excluding a stage that is included as a macro will exclude that stage and not call it's dependents, so the following calls are equivalent::
			
			$ SEDMsim -stage *stage
			$ SEDMsim *none
			
		
	The commonly used stages in :program:`SEDMsim` are
	
	=====================  ========================================
	  Stage                Description                           
	=====================  ========================================
	  ``*all``              Run all default stages                 
	  ``*cached-only``      Use cached subimages to construct final image                  
	  ``*none``             Run no stages                          
	  ``*setup-blank``		Set up a blank image.
	  ``*crop``             Crop Final Image                       
	  ``*add-noise``        Add Dark/Bias noise to image           
	  ``*transpose``        Transpose the image                    
	  ``*save``             Save image to file                     
	  ``*simple-source``    Use a simple, centered source object   
	  ``*sky-source``       Use only sky spectrum                  
	  ``*flat-source``      Use a constant value source            
	  ``*line-source``      Use a calibration lamp source          
	  ``*setup``            System Setup                           
	  ``*plot``             Run all stages that make plots         
	=====================  ========================================
	
	Stages are registered by the :meth:`AstroObject.AstroSimulator.Simulator.registerStage` method.
	
 .. option:: [configurations]
	
	Configuration options override defaults set up in :class:`SEDMachine.Simulator.SEDSimulator`. As such, they are useful quick changes to a configuration.
	
	
	 ===================== ============================================
	   Options               Description
	 ===================== ============================================
	   ``-d``                enable debugging messages and plots
	   ``-D``                Debug, Limit lenslets (50,start from 2100)
	   ``-S``                Debug, Limit lenslets (5,start from 2100)
	   ``-T``                Limit lenslets (50,start from 2100)
	   ``-M``                Limit lenslets (500,start from 1000)
	   ``-N``                Limit lenslets (500,start from 1000)
	   ``-A``                Debug, Single lenslets (start from 2100)
	 ===================== ============================================
	
	.. Note::
		This means if you dump the configuration file (see :option:`--dump`) and use that directly as your new configuration file, these options will have no effect. Therefore, it is advisable that your configuration file contain the minimum amount of detail to override the default values set in the program. However, if you wish to use these options, but always disable debug, for example, you could disable debug in your configuration file. This will make none of these flags enable debugging.
		
	
 .. option:: -h, --help
	
	Get help about the :program:`SEDMsim` command. Help will list commonly used stages, optional arguments, and configuration items.
	
	.. Note ::
		To get a full list of stages available, use :option:`--stages`
		
 .. option:: --version
	
	Print the program version.
	
 .. option:: --cf file.yaml
	
	Set the name of the configuration file to use. By default, the configuation file is called `SED.main.config.yaml`. This will be the filename used for dumping configurations with the :option:`--dump` command (Dump will append the extension ``-dump.yaml`` onto the filename to prevent overwriting exisitng configurations)
	
 .. option:: --dry-run
	
	Traverse all of the stages to be run, printing them as the program goes, but do not run any stages.
	
 .. option:: --stages
	
	Print all stages registered in the simulator. Any stage listed in the output of this function can be run.
	
 .. option:: --dump
	
	Dump the configuration to a file. Filenames are the configuration file name with ``-dump.yaml`` appended.

.. _Stages:

:program:`SEDMsim` Stages
-------------------------

The following are the most commonly used :program:`SEDMsim` stages in a production environment.

 .. describe:: *all
	
	Complete a full simulation.
	
 .. describe:: *cached-only
    
	Complete only the parts of the simulation conducted after the ``*place`` stage. I.e. assemble a full image from pre-rendered spectra stored in the ``Caches/`` folder.

 .. describe:: *none
	
	Do no stages.
	
 .. describe:: *setup-blank
 	
	Make a basic, empty image.
	
 .. describe:: *crop
 	
	Crop the image to the appropriate CCD size.
 
 .. describe:: *add-noise
 
 	Add read noise and dark current appropriate for the image.
	
 .. describe:: *transpose
 
    Account for the difference between x,y ordered arrays and c-ordered (y,x) arrays that are used in Numpy.
	
 .. describe:: *save
 
    Save the image to a ``.fits`` file.

 .. describe:: *simple-source
	
	Create a simple source for use in the simulation. The simple source geometry is hardcoded into the program. To examine the geometry of this source, use::
		
		$ SEDMsim *simple-source *plot-p-geometry
		
	
	See :meth:`SEDMachine.Simulator.SEDSimulator.setup_simple_source`
	
 .. describe:: *flat-source
 
 	Use only a flat, uniformly applied source. Source will use the system throughput and quantum efficiency. To disable system throughput interpretation, use ``-apply-qe``.
	
 .. describe:: *line-source
 
 	Use a list of lines as the source. The list of lines should be at minimum a list of wavelengths in meters. If desired, the list of lines can be two columns, the first for the position of the line in meters, the second for the height of the peak of the line in Joules per meter (similar to ergs per angstrom), i.e. F-Lambda units.
	
 .. describe:: *sky-source
 
 	Use only the sky, uniformly applied, in the image.

 .. describe:: *setup
 	
	Does all of the required methods to ready the simulator for action. This stage is a good test of the configuration and basic opertaion of the system without much computational expense.

 .. describe:: *plot
	
	Make plots about the system. Plots made are:
	
	- The pupil x,y positions of lenslets
	- The camera x,y positions of lenslets
	- Arc-distance (delta) in px vs. x in pixels. ``%(Partials)s/Instrument-%(num)04d-Delta-Distances%(ext)s``
	- Wavelength vs. effective sampling resolution. ``%(Partials)s/Instrument-%(num)04d-WL%(ext)s``
	

.. _SEDM_Configuration:

:program:`SEDMsim` Configuration Files
--------------------------------------

:program:`SEDMsim` configuration files are YAML files which contain a dictionary structure. All values in the YAML files are basic yaml, and contain no python-specific directives. To find out what the default or current configuration is, use the :option:`--dump` command. The file produced from this will contain a YAML structure for the configuration in use when the system started up. The various directives in the configuration file are described below.

In general, the configuration is minimally set. If you don't know what a value does, don't set it! The defaults should do just fine. If you need to see what a minimal configuration looks like, examine the :program:`SEDMsetup` ``*config-file`` stage.

The configuration is documented with Dot-syntax. Each dot represnets a new sub-level. As such, **Instrument.Scatter.Kernel.[].mag** would be written as
::
	
	Instrument:
	  Scatter:
	    Kernel:
	      A:
	        mag: value

where ``[]`` is used as a substitue for some key, ``A`` in this case.


General Configuration
~~~~~~~~~~~~~~~~~~~~~

 .. object:: Configurations
	
	Files used to configure the system.
	
 .. object:: Configurations.Main
		
	The name of the configuration file.
		
 .. object:: Configurations.This
		
	The most recently selected and read configuration file, and the one that will be loaded next.
		
 .. object:: Debug
	
	If true, will halt the simulation when it comes across exceptions. False will run the simulator, trying to recover from as many exceptions as possible.
	
 .. object:: Dirs
 	
 	Directories used in the simulation environment. This configuration is used to set up the environment in :program:`SEDMsetup` and then used to find data files, and store other files during the simulation.
	
 .. object:: Dirs.Caches
	
	A directory to store files which will be used internally by the simulator. These files may provide useful information about the work of the simulator and are generally written in standard formats. This directory will be read when using the ``*cached-only`` stage. Files in here should generally only be created by the :program:`SEDMsim` program, to ensure a consistent and expected format. The contents of this directory can be deleted at any point, and with the exception of the special command ``*cached-only``, its contents will be regenerated on the next run.
	
 .. object:: Dirs.Output
	
	The output directory for the system.
	
 .. object:: Dirs.Logs
	
	The logging directory for the system. Logs are stored here by default, and rotated at midnight each night. Logging can be quickly disabled by removing this folder.
	
 .. object:: Dirs.Partials
		
	Files generated during the simulation which are useful introspections of the process, but not necessary outputs. Plots, and readouts of data are stored here.
		
	

Instrument Configuration
~~~~~~~~~~~~~~~~~~~~~~~~

 .. object:: Instrument
	
	Configuration related to the characteristics of the instrument iteslf.
	
 .. object:: Instrument.Thpt
	
	Data about the throughput of the system.
	
 .. object:: Instrument.Thpt.File
	
	The numpy throughput calculation file. This file is provided as-is, and no customization of the file is available right now.
	
 .. object:: Instrument.Thpt.Type
	
	The type of throughput system to use. The value ``prism_pi`` uses the characteristics of the three-glass prism in SEDMachine and the PI CCD.
	
 .. object:: Instrument.Cameras
	
	A dictionary of camera settings
	
 .. object:: Instrument.Cameras.[]
	
	The identifying key for each camera type.
	
 .. object:: Instrument.Cameras.[].DC
	
	Camera dark current, per exposure second, in counts.
	
 .. object:: Instrument.Cameras.[].RN
	
	Camera read noise, per exposure, in counts.
	
 .. object:: Instrument.Cameras.[].readtime
	
	Camera read time in seconds
	
 .. object:: Instrument.Cameras.Selected
	
	The name of the selected camera.
	
 .. object:: Instrument.Lenslets
	
	Information about the lenslet array.
	
 .. object:: Instrument.Lenslets.radius
	
	Raidus of each lenslet in mm
	
 .. object:: Instrument.Lenslets.rotation
	
	Rotation of the lenslet, relative to CCD columns, in degrees.
	
 .. object:: Instrument.Lenslets.strict
	
	Determines whether lenslets are included in a strict fashion, or whether any possible lenslet should be calculated. Should only be turned off during advanced debugging, as in general, the strict validity checks will not eliminate any lenslet which has any effect on the system. See :meth:`SEDMachine.Objects.Lenslet.valid`.
	
 .. object:: Instrument.PSF
	
	Data about the instrument PSF
	
 .. object:: Instrument.PSF.size
	
	The size of the PSF in pixels.
	
 .. object:: Instrument.PSF.stdev
	
	The standard deviation of the PSF (if a gaussian is used.)
	
 .. object:: Instrument.Scatter
	
	Settings to do with scattered light generation.
	
 .. object:: Instrument.Scatter.Amplifier
	
	A scalar gain setting for scattered light.
	
 .. object:: Instrument.Scatter.FFT
	
	Boolean, whether the system should use a fast-fourier-transform algorithm.
	
 .. object: Instrument.Scatter.Kernels
	
	A list of kernels for the scattered light profile. Each kernel should have the following:
	
 .. object:: Instrument.Scatter.Kernel.[].mag
	
	Maximum height of the kernel.
	
 .. object:: Instrument.Scatter.Kernel.[].stdev
	
	Standard deviation of the kernel.
	
 .. object:: Instrument.Scatter.Kernel.[].type
	
	"Gaussian" for gaussian kernels (the only ones supported right now.)
	
 .. object:: Instrument.Tel
	
	Information about the telescope.
	
 .. object:: Instrument.Tel.area
	
	Telescope mirror area in square mm.
	
 .. object:: Instrument.Tel.dispfitorder
	
	Fit order used when fitting ellipse size from the ray-trace for each lenslet.
	
 .. object:: Instrument.Tel.ellipse
	
	Whether elliptical shape data is included in the ray-trace file.
	
 .. object:: Instrument.Tel.obsc.{ratio,px}
	
	Size of the telescope's central obscuration. If elliptical images are used, should be specified as a fraction of the radius of the telescope image/primary mirror diameter, using the ``ratio`` keyword. If not, should be specified in ``mm`` or ``px`` keywords in absolute terms.
	
 .. object:: Instrument.Tel.raidus.px
	
	Radius of the primary mirror image in pixels.
	
 .. object:: Instrument.ccd
	
	Information about the CCD
	
 .. object:: Instrument.ccd.size.px
	
	Size of the ccd in px.
	
 .. object:: Instrument.convert.pxtomm
	
	The number of pixels per mm.
	
 .. object:: Instrument.density
	
	The amount of oversampling used to generate spectra.
	
 .. object:: Instrument.dispfitorder
	
	Wavelength-to-arclength function fit order.
	
 .. object:: Instrument.eADU
	
	Electrons per photon.
	
 .. object:: Instrument.files
	
	Files used in the simulator.
	
 .. object:: Instrument.files.encircledenergy
 	
	An encircled energy function, showing percentage of energy encircled and distance from central pixel.
	
 .. object:: Instrument.files.lenslets
 	
	Ray trace file for the system. See :meth:`SEDMachine.Simulator.SEDSimulator.setup_lenslets`
	
 .. object:: Instrument.image
 	
	Information about the image plane. The image plane needs to contain all photons from all traces to prevent boundary errors.
	
 .. object:: Instrument.image.pad.mm
 	
	Padding between ray traces and the edge of the image plane.
	
 .. object:: Instrument.image.size.mm
	
	Size of the image plane
	
 .. object:: Instrument.padding
 	
	Pixels (full-size) of padding around each spectrum for individual specturm generation. Must be **larger** than the convolved PSF and telescope image in order to capture all photons in each sub-image.
	
 .. object:: Instrument.wavelengths
 	
	Wavelenght settings for arbitrary spectra (not final image) generation used for debugging.
	
 .. object:: Instrument.wavelengths.max
	
	Longest wavelegnth included.
	
 .. object:: Instrument.wavelengths.min
	
	Shortest wavelength sampled.
	
 .. object:: Instrument.wavelengths.resolution
 	
	Constant sampling resolution.
	

Observation and Source Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Logging and other Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 .. object:: Caches
	
	Cache object Filenames
	
 .. object:: Caches.CONV
		
	Cache filename for the convolution of the telescope image and the PSF.
		
 .. object:: Caches.PSF
		
	Cache filename for the instrument PSF
		
 .. object:: Caches.Telescope
		
	Cache filename for the telescope image
		
 .. object:: Caches.config
		
	Cache filename for the configuration
		
 .. object:: Caches.const
		
	Cache filename for the set of constants used in the program.
		
	
