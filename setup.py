from setuptools import setup, find_packages
setup(
    name = "SEDMachineSimulator",
    version = "0.1.3p2",
    packages = find_packages(
        exclude=['Tests.*'],
        ),
    include_package_data = True,
    exclude_package_data = {'': ['.gitignore','bump-version.sh','distribute.sh'], 'Docs/build':['*']},
    install_requires = ['pyfits>=2.4','numpy>=1.6','scipy>=0.9','matplotlib>=1.0','AstroObject>=0.2,<0.3.0a0'],
    dependency_links = ['https://github.com/alexrudy/AstroObject/zipball/v0.2.8#egg=AstroObject-0.2.8'],
    test_suite = 'Tests',
    entry_points = {
        'console_scripts' : ['SEDMsim = SED.Simulator:run',],
    },
)
