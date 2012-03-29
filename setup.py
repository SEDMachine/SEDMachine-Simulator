from distribute_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages
setup(
    name = "SEDMachineSimulator",
    packages = ['SEDMachine','SEDTools'],
    package_data = {'':['VERSION','README.md','LICENCE','*.yaml'],'SEDTools':['Data/*.dat','Data/*.npy']},
    version = "0.3.2",
    install_requires = ['pyfits>=2.4','numpy>=1.5','scipy>=0.9','matplotlib>=1.0','AstroObject>=0.3.3'],
    dependency_links = ['https://github.com/alexrudy/AstroObject/zipball/v0.3.3#egg=AstroObject-0.3.3'],
    test_suite = 'Tests',
    entry_points = {
        'console_scripts' : ['SEDMsim = SEDMachine.Simulator:run', 'SEDMsetup = SEDTools.setup:run'],
    },
)
