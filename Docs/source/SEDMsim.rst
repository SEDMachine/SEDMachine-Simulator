.. program:: SEDMsim
.. _SEDMsim:

:program:`SEDMsim` Program for running simulations
============================================================

The master simulator program is a command-line interface to the :meth:`AstroObject.AstroSimulator.Simulator.run` method. Below are the major command line components.


.. option:: {stages}
	
	The stages option specifies individual stages for the program to run. You must specify at least one stage to run in the simulator. By default, two basic stages are provided, ``*all`` and ``*none``. The default simulation is performed by `*all`. To test the simulator without running any stages (for example, to test :meth:`AstroObject.AstroSimulator.Simulator.registerFunction` functionality), us the `*none` stage to opertate without using any stages.
	
	Stages are called with either a ``*``, ``+`` or ``-`` character at the beginning. Their resepctive actions are shown below.
	
	========= ============ ================================
	Character  Action      Description
	========= ============ ================================
	``*``     Include      To include a stage, use ``*stage``. This will also run the dependents for that stage.
	``-``     Exclude      To exclude a stage, use ``-stage``. This stage (and it's dependents) will be skipped.
	``+``     Include-only To include a stage, but not the dependents of that stage, use ``+stage``.
	========= ============ ================================
	
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