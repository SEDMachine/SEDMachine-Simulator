.. program:: SEDMsetup
.. _SEDMsetup:

:program:`SEDMsetup` Program managing environments
==================================================

This program is a command-line interface to set-up a simulation environment for the :program:`SEDMsim` tool.

:program:`SEDMsetup` Usage
--------------------------

Usage Statement::
	
	SEDMsetup [ options ] { stages }
	
.. option:: {stages}
	
	The stages option specifies individual stages for the program to run. You must specify at least one stage to run in the simulator.
	The commonly used stages in :program:`SEDMsetup` are
	
	=====================  ========================================
	  Stage                Description                           
	=====================  ========================================
	  ``*all``              Run all default stages                 
	  ``*none``             Run no stages                          
	  ``*config-file``     Generate a basic configuration file
	=====================  ========================================
	
	Stages are registered by the :meth:`AstroObject.AstroSimulator.Simulator.registerStage` method.

.. option:: -h, --help
	
	Get help about the :program:`SEDMsetup` command. Help will list commonly used stages, optional arguments, and configuration items.
	
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
	

:program:`SEDMsetup` stages
---------------------------

The basic, important stages are:

 .. object:: *all
 	
	Produces directories for the system. Loads spectral data files and basic ray trace files which are packaged with the simulator.
	
 .. object:: *config-file
 	
	Generates a simple configuration file which correctly uses the packaged data. The basic configuration file contains all instrument-specific values (as opposed to algorithm or operationally specific values). To get a full configuration file (which includes all values...) use::
	
	$ SEDMsim --dump *none
	

:program:`SEDMsetup` configuration
----------------------------------

This program respects the basic configuration items of :program:`SEDMsim`