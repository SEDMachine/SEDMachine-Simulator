Simulation Runner
=================

The command line runner is an object for handling command line interaction, and which calls and controls all of the various command line interaction procedures. To use the command line runner in your own script, use ::
    
    if __name__ == '__main__':
        Sim = Simulator()
        Sim.run()
    

This will set up the command line runner, and activate it. In general, the command line simulator runner is not pluggable. The architecture has not been designed to easily include other simulation components. As such, if you wish to expand on the command line runner, it will take a bit of hacking. The functions in the command line runner are documented below. In general, only :func:`Simulator.run()` is designed to be a publicly accessible function.

.. Note :: Version 0.2 will use :mod:`AstroObject` v0.3.0, which contains a simulator model that is highly pluggable. Plugging the simulator will make more sense then.

.. autoclass::
    SED.Simulator
    :members: