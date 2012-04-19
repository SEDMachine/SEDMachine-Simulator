from distribute_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

from distutils.command.build_py import build_py as du_build_py
from distutils.core import Command

from SEDMachine.version import version as versionstr

class Version(Command):
    description = "Print the module version"
    user_options = []
    def initialize_options(self):
        pass
        
    def finalize_options (self):
        pass
        
    def run(self):
        print 'version',versionstr
        return
        

 

#custom build_py overwrites version.py with a version overwriting the revno-generating version.py
class SED_build_py(du_build_py):
    def run(self):
        from os import path
        res = du_build_py.run(self)
        
        versfile = path.join(self.build_lib,'SEDMachine','version.py')
        print 'freezing version number to',versfile
        with open(versfile,'w') as f: #this overwrites the actual version.py
            f.write(self.get_version_py())
        
        return res
        
    def get_version_py(self):
        import datetime
        from SEDMachine.version import _frozen_version_py_template
        from SEDMachine.version import version,major,minor,bugfix,dev
        
        
        timestamp = str(datetime.datetime.now())
        t = (timestamp,version,major,minor,bugfix,dev)
        return _frozen_version_py_template%t
        

SEDpkgs = find_packages(exclude=['Tests'])
AstroObjectReq = "0.5-a1"
AstroObjectDep = "AstroObject>=" + AstroObjectReq
AstroObjectVer = "0.5-a1"
AstroObjectURL = "https://github.com/alexrudy/AstroObject/zipball/v%(ver)s#egg=AstroObject-%(ver)s" % { 'ver' : AstroObjectVer}

class AstroObjectSourceURL(Command):
    description = "Print the AstroObject URL"
    user_options = []
    def initialize_options(self):
        pass
        
    def finalize_options (self):
        pass
        
    def run(self):
        print 'AstroObject:',AstroObjectURL
        return


setup(
    name = "SEDMachineSimulator",
    packages = SEDpkgs,
    package_data = {'':['VERSION','README.md','LICENCE','*.yaml'],'SEDTools':['Data/*.dat','Data/*.npy','*.yaml']},
    version = versionstr,
    install_requires = ['pyfits>=2.4','numpy>=1.5','scipy>=0.9','matplotlib>=1.0','shapely>=1.2.14',"PyYAML>=3.10",AstroObjectDep],
    dependency_links = [AstroObjectURL],
    test_suite = 'Tests',
    author = "Alexander Rudy",
    author_email = "dev@alexrudy.org",
    entry_points = {
        'console_scripts' : ['SEDMsim = SEDMachine.Simulator:run', 'SEDMsetup = SEDTools.setup:run','RCpipeline = RCPipeline.pipeline:main','RCsetup = RCPipeline.tools:main'],
        'distutils.commands' : ['version = SEDMachine.version:show',],
    },
    cmdclass = {
        'build_py' : SED_build_py,
        'version' : Version,
        'aosource': AstroObjectSourceURL,
    },
)
