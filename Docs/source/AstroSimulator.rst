.. module:: AstroObject.AstroSimulator

Simulator and Complex Task Management :mod:`AstroSimulator`
***********************************************************

The Simulator module is designed to provide a high level, command-line useful interface to large computational tasks. As the name suggests, Simulators often do a lot of programming work, and do so across many distinct "stages", whcih can be configured in any way the user desires. Each stage should operate as an independent, no argument function. Stages may rely on a variety of dependencies to ensure that other actions happen first. As well, stages may substitute for others to provide alternative functionality. As stages generally run in their appropriate order, and then follow dependency tress downwards, ensuring that the simulator does a minimum amount of work for each task is relatively easy.

The simulator API shown below is for reference when understanding and expanding the :ref:`SEDMsim` tool.

.. autoclass::
    AstroObject.AstroSimulator.Simulator
    
	.. automethod:: registerStage
	
	.. automethod:: registerConfigOpts
	
	.. automethod:: registerFunction
	
	.. automethod:: run
	
	.. automethod:: startup
	
	.. automethod:: do
	
	.. automethod:: execute
