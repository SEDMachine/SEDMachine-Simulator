from distribute_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages
setup(
    name = "SEDMachineSimulator",
    packages = ['SEDMachine','SEDSpec'],
    package_data = {'':['VERSION','README.md','LICENCE'],'SEDspec':['*.npy','*.dat']},
    version = "0.2.0a1",
    install_requires = ['pyfits>=2.4','numpy>=1.5','scipy>=0.9','matplotlib>=1.0','AstroObject>=0.3.0a2'],
    dependency_links = ['https://github.com/alexrudy/AstroObject/zipball/v0.3.0a1#egg=AstroObject-0.3.0a1'],
    test_suite = 'Tests',
    entry_points = {
        'console_scripts' : ['SEDMsim = SEDMachine.Simulator:run',],
    },
)
