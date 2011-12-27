Command Line Interface
======================

The SED Machine Data Simulator is designed to be run as a command line tool. It has a variety of configuration files, as well as some command line options, which can be used to control the Simulator's behavior. This section will document their use. The underlying python code is also documented, and can be viewed in the following sections.

The command line tool is by default called ``SED.py``. As such, it is normally run as ``./SED.py``.

For basic help, run::
    
    $ ./SED.py -h
    usage: ./SED.py [-D | -T | -E | -F] [arguments] subcommand
    
    Utilities for simulating the SEDMachine
    
    optional arguments:
      -h, --help            show this help message and exit
      -f label              label for output image
      --version             show program's version number and exit
      -D, --dev             equivalent to --debug --easy
      -T, --test            equivalent to --no-cache --easy
      -E, --easy            uses simple configuration for testing the simulator
      -F, --full            disables --debug and forces normal deployment
                            operation
      --plot                Enable debugging plots
      --no-cache            ignore cached data from the simulator
      -n N                  limit the run to N lenslets
      -o I                  start filling from lenslet in array position I
      -d, --debug           enable debugging messages and plots
      --config file.yaml    use the specified configuration file
      
    sub-commands:
      {all,poscache,startup,postest,source}
          all                 run the full simulator
          startup             run the initalization and caching for the system
          source              run the source creation and resolution routines
        
    This is the Command Line Interface to the SEDMachine Simulator Program.
    

The ``SED.py`` script is generally structured with global options first, then a subcommand which controls the specific behavior of the system, and finally, arguments which pertain specifically to that subcommand.

The Global Options and each subcommand is documented in turn below.

Mode Options: ``-D``, ``-T``, ``-E``, ``-F``
--------------------------------------------
This set of options declares simple configuration settings for easy use of the program

:option:`-F`

By default, the ``-F`` flag is set, which provides no simple configuration options and will not override configuration files in any way. 

:option:`-E`

The ``-E`` option is for Easy Mode, which uses a uniform blackbody source at 5000K to illuminate the entire lenslet array, and then selects only 5 lenslets for dispersion and placement. 

:option:`-D`:option:`-T`

Both the ``-D`` and ``-T`` options imply the ``-E`` option settings, but the ``-D`` option (for development) enables debugging messages in logs, and the ``-T`` option disables caching (so that you can ensure the script will run correctly on a fresh system).

The ``all`` Subcommand
----------------------

The ``SED.py all`` command is designed to run a full simulation. In this mode, the system will start up, initialize all of the required components, and produce a fully simulated image. In order to operate in this mode, the instrument and source must both be correctly configured. This subcommand can take an optional ``-c`` parameter which can be used to override the configuration file used for the source.

::
    
    $ ./SED.py all -h
    usage: ./SED.py [-D | -T | -E | -F] [arguments] all [-c config ]
    
    Runs the simulator with source configuration provided
    
    optional arguments:
      -h, --help  show this help message and exit
      -c config
    

.. Note:: Currently, the system only supports uniform sources. See the Source Configuration File


.. Note:: If the program is ever updated, this mode will be updated to include the latest simulation components.


The ``startup`` Subcommand
--------------------------

The ``SED.py startup`` command is designed to simply test the initialization of all features in the Instrument simulator (but not the source or noise simulation). It is useful to test instrument configurations for validity, as an invalid configuration should strike an error before completing the startup process.

The ``source`` Subcommand
-------------------------

The ``SED.py source`` command runs through the initialization of all features in the Instrument and Source simulators, as well as generating the noise masks. As such, running this command is a good indicator that all of the configuration validates and can run.

Other Command Options
---------------------

The command has a few remaining command line switches which can be used to over-ride configuration settings at runtime. The options are discussed below.

:option:`-f`

The label for the output file. This label will be used in constructing the full filename, which will also include a date.

:option:`--plot`

Enables diagnostic plotting. Plotting makes the script much slower (it uses LaTeX enabled :mod:`MatPlotLib`), but will produce many useful plots if :option:`--debug` is enabled.

:option:`--no-cache`

Disables any caching functions, forcing the script to generate all information newly on this run.

:option:`-d --debug`

Enables debugging mode

:option:`-n N`

Limits the number of spectra to be placed in the final image to ``N``.

:option:`-o I`

Sets the first spectrum to be included to ``I``

:option:`--config`

Sets the configuration file to the given filename.

:option:`--dump-config`

Writes all of the configuration files out to new files. This is a good way to generate template configuration files.

Configuration Files
===================

All of the configuration files are YAML, and are used to override the default configuration values. As such, configuration files need only include those values which differ from the defaults, and do not have to include all possible configuration terms. 

A very short configuration file could be::
    
    Cache: true
    
This valid configuration file would use the default value for all variables, except caching, which would be enabled.

To generate example configuration files (with all of the default variables) for use as templates, use the :option:`--dump-config` option::
    
    $ SED.py --dump-config all
    
This will generate files which end in ``.dump.yaml`` which are the dumped configuration files from the system.

A basic script configuration looks like::
    
    Cache: true
    CacheFiles:
      Instrument: SED.instrument
    Configs:
      Instrument: SED.instrument.config.yaml
      Source: SED.source.config.yaml
      This: SED.script.config.yaml
    Debug: false
    Dirs:
      Caches: Caches/
      Images: Images/
      Logs: Logs/
      Partials: Partials/
    Lenslets:
      number: 5
      start: 2150
    Plot: false
    

Script Configuration ``SED.script.config.yaml``
-----------------------------------------------

This file configures the operation of the simulator script, and its general behavior. This file should not include any scientific specifics.

::
    
    Cache: true
    

This key enables or disables caching. Caching is used to store telescope images and wavelength positions. In the future, it may store source information as well. The caching engine tries to be intelligent, and notice altered configurations, so it is generally safe to enable caching. However, if there are changes to the lenslet specification, or dispersion specification, the caching engine will not notice (as these changes are hard to detect in outcome variables) and so when these files change, caching should be disabled.

::
    
    CacheFiles:
      Instrument: SED.instrument
    

This directive sets the base filename for caching files. By default, caching files can be found in ``Caches/`` and each file name will start with this string. It is relatively arbitary, and used only internally. However, if you wish to keep two different systems separate, you could change the name of the caching files for each system. As such, you would not have to regenerate caches every time you switched systems.

::
    
    Configs:
      Instrument: SED.instrument.config.yaml
      Source: SED.source.config.yaml
      This: SED.script.config.yaml
    

This directive sets the names of all of the configuration files. The ``This`` parameter is self-referntial, i.e. refers to the file name of the script configuration file.

::
    
    Debug: false
    

This directive enables or disables debugging output. Debugging must be enabled for plotting output.

::
    
    Dirs:
      Caches: Caches/
      Images: Images/
      Logs: Logs/
      Partials: Partials/
    

These are the directories used by the script to store various materials. The ``Caches`` directory holds cached data. This data can be deleted at any time, and will be regenerated if necessary. The ``Images`` directory holds output images (i.e. final, simulated CCD frames). The ``Paritals`` directory holds partially generated items which are usually only populated with debugging enabled. This directory is used for text-file readouts of data, as well as intermediate plots and other diagnostic information. The ``Logs`` directory simply holds log files generated by the system, which often contain useful diagnostic data about the simulation, and can be consulted if the simulation crashes.

::
    
    Lenslets:
      number: 5
      start: 2150
    

The lenslet directives are useful for limiting the number of lenslets used for placing spectra. The ``number`` is simply the number to place, and the ``start`` is an index offset from which to start, allowing you to only place five lenslets in the middle of the ccd (this example) for testing purposes.

::
    
    Output:
      Format: fits
      Label: Generated
    

The output directive controls information about the final generated image. ``Label`` is the text to be included in the image name, and ``Format`` is the extension to be used.

.. Note:: Currently, the program only generates FITS format files, so fits and fit are the only format parameters that make sense.

::
    
    Plot: false
    

This parameter controls debugging plots, which make the system a lot slower but provide useful intermediate step information.

::
    
    Source:
      Temp: 5000
      Type: BlackBody
    

The source directive provides a way to set systemwide defaults for sources used in simulations. These values will be used in place of the program defaults for the source system, but will be overridden by a specific source configuration file.

Instrument Configuration ``SED.instrument.config.yaml``
-------------------------------------------------------
This file is used to configure the instrument model.

::
    
    convert:
      pxtomm: 0.0135
      

Units are handled intelligently in configuration files. This value sets the conversion between pixels and mm. Any value that is a distance can then be provided with either the ``px`` keyword or the ``mm`` keyword. The program will automatically calculate the conversion between the two. To disable conversion calculation for a specific value, add a ``calc: false``

::
    
    bias: 20
    

This parameter sets the bias level in mean counts per pixel. Bias is then generated using poisson statistical random noise.

::
    
    ccd_size:
      px: 2048
    

    dark: 20
    density: 5
    exposure: 120
    files:
      dispersion: Data/dispersion_12-10-2011.txt
      encircledenergy: Data/encircled_energy_4nov11.TXT
      lenslets: Data/xy_17nov2011_v57.TXT
    gain: 1.0e-06
    image_size:
      mm: 40.0
    logging:
      console:
        enable: false
        format: '... ...%(message)s'
        level: false
      file:
        enable: true
        filename: SED-2011-12-03
        format: '%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s'
    padding: 5
    plot_format: .pdf
    psf_size:
      px: 0
    psf_stdev:
      px: 1.0
    tel_obsc:
      px: 0.2
    tel_radii:
      px: 1.2
    

Source Configuration ``SED.source.config.yaml``
-----------------------------------------------

Logging Configuration
---------------------
The logging directives can be placed in any configuration file, and will configure logging for that item.

The logging configuration directives are::
    
    logging:
      console:
        enable: false
        format: '... ...%(message)s'
        level: false
      file:
        enable: true
        filename: SEDSource-2011-12-03
        format: '%(asctime)s : %(levelname)-8s : %(funcName)-20s : %(message)s'
    

.. Note:: These directives will not be well documented here, but can be adjusted if you need. In general, the ``enable`` directive should be all you really need to enable and disable various forms of logging.




