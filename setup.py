from distribute_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

from distutils.command.build_py import build_py as du_build_py
from AstroObject.version import version as versionstr

#custom build_py overwrites version.py with a version overwriting the revno-generating version.py
class ao_build_py(du_build_py):
    def run(self):
        from os import path
        res = du_build_py.run(self)
        
        versfile = path.join(self.build_lib,'AstroObject','version.py')
        print 'freezing version number to',versfile
        with open(versfile,'w') as f: #this overwrites the actual version.py
            f.write(self.get_version_py())
        
        return res
        
    def get_version_py(self):
        import datetime
        from AstroObject.version import _frozen_version_py_template
        from AstroObject.version import version,major,minor,bugfix,dev
        
        
        timestamp = str(datetime.datetime.now())
        t = (timestamp,version,major,minor,bugfix,dev)
        return _frozen_version_py_template%t
        


setup(
    name = "SEDMachineSimulator",
    packages = ['SEDMachine','SEDTools','SEDTools.Data'],
    package_data = {'':['VERSION','README.md','LICENCE','*.yaml'],'SEDTools':['Data/*.dat','Data/*.npy','*.yaml']},
    version = versionstr,
    requires = ['pyfits>=2.4','numpy>=1.5','scipy>=0.9','matplotlib>=1.0','AstroObject>=0.3.4'],
    dependency_links = ['https://github.com/alexrudy/AstroObject/zipball/v0.3.4#egg=AstroObject-0.3.4'],
    test_suite = 'Tests',
    author = "Alexander Rudy",
    author_email = "dev@alexrudy.org",
    entry_points = {
        'console_scripts' : ['SEDMsim = SEDMachine.Simulator:run', 'SEDMsetup = SEDTools.setup:run'],
    },
    cmdclass = {
        'build_py' : ao_build_py,
    }
)
