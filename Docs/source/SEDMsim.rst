.. program:: SEDMsim
.. _SEDMsim:

:program:`SEDMsim` Program for running simulations
==================================================

The master simulator program is a command-line interface to the :meth:`AstroObject.AstroSimulator.Simulator.run` method. Below are the major command line components.

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
	  ``*none``             Run no stages                          
	  ``*plot``             Run all stages that make plots         
	  ``*simple-source``    Use a simple, centered source object   
	  ``*setup``            System Setup                           
	  ``*sky-source``       Use only sky spectrum                  
	  ``*flat-source``      Use a constant value source            
	  ``*line-source``      Use a calibration lamp source          
	  ``*crop``             Crop Final Image                       
	  ``*add-noise``        Add Dark/Bias noise to image           
	  ``*transpose``        Transpose the image                    
	  ``*save``             Save image to file                     
	  ``*cached-only``      Use cached subimages to construct final image                  
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

.. _Configuration:

:program:`SEDMsim` Configuration Files
======================================

:program:`SEDMsim` configuration files are YAML files which contain a dictionary structure. All values in the YAML files are basic yaml, and contain no python-specific directives. To find out what the default or current configuration is, use the :option:`--dump` command. The file produced from this will contain a YAML structure for the configuration in use when the system started up. The various directives in the configuration file are described below.

- Caches:
	- CONV: SED.conv.npy
	- PSF: SED.psf.npy
	- Telescope: SED.tel.npy
	- config: SED.config.yaml
	- const: SED.const.yaml
- Configurations:
	- Main: SED.main.config.yaml
	- This: SED.main.config.yaml
- Debug: false
- Dirs:
	- Caches: Caches
	- Images: Images
	- Logs: Logs
	- Partials: Partials
- Instrument:
	- Thpt:
	    - File: SEDSpec2/Data/thpt.npy
	    - Type: prism_pi
	- bias: 20
	- camera: PI
	- ccd_size:
	    - px: 2048
	- convert:
	    - pxtomm: 0.0135
	- density: 5
	- dispfitorder: 2
	- eADU: 0.03802
	- files:
		- dispersion: Data/dispersion_12-10-2011.txt
		- encircledenergy: Data/encircled_energy_4nov11.TXT
		- lenslets: Data/xy_17nov2011_v57.TXT
	- gain: 5
	- image_pad:
		- mm: 0.1
	- image_size:
		- mm: 40.0
	- lenslets:
		- radius: 0.00245
		- rotation: 27.0
	- padding: 5
	- plot: false
	- psf_size:
	    - px: 2.4
	- psf_stdev:
		- px: 1.0
	- scatter:
		- radius: 800
		- wavelength: 4.5e-07
	- tel_area: 16417.8
	- tel_obsc:
	    - px: 0.2
	    - ratio: 0.1
	- tel_radii:
	    - px: 1.2
	- wavelengths:
	    - max: 9.3e-07
	    - min: 3.7e-07
	    - resolution: 100
- Lenslets: {}
- Observation:
	- Moon:
	    - Phase: 6
	- Sky:
		- Atmosphere: Atmosph
		- Files:
			- HansuchikUVES: SEDSpec2/HansuchikUVES.fits
			- Massey: SEDSpec2/MasseySky.fits
			- PALext: SEDSpec2/atmosphere.fits
			- Quimby: SEDSpec2/Quimby.fits
			- TurnroseSKY: SEDSpec2/TurnroseSKY.fits
		- Use: TurnroseSKY
	- airmass: 1
	- exposure: 1200
	- number: 3
- Output:
	- Format: fits
	- Label: Generated
- Plots:
	- format: .pdf
- Source:
	- CubeName: Data/CUBE.fits
	- Filename: Data/SNIa.R1000.dat
	- Flat:
		- value: 1.0e-06
	- LineList: null
	- PXSize:
		- mm: 0.005
	- PreAmp: 0.001
	- Rotation: 0.7853981633974483
	- Sample_Lenslet: 2000
	- WLCal:
		- List: Data/Lines.dat
	    - sigma: 1.0e-09
	    - value: 100000000.0
- logging:
	- console:
		- enable: true
	    - format: '%(levelname)-8s... %(message)s'
	    - level: 20
	    - message: '...%(message)s'
	- file:
	    - dateformat: '%Y-%m-%d-%H:%M:%S'
	    - enable: true
	    - filename: SEDMachine
	    - format: '%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s'
	    - level: null
	- growl:
	    - enable: false
	    - level: null
	    - name: SED Machine Simulator
