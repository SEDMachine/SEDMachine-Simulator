from setuptools import setup, find_packages
setup(
    name = "SEDMachineSimulator",
    version = "0.1.2",
    packages = find_packages(
        exclude=['Tests.*'],
        ),
    exclude_package_data = {'': ['.gitignore','bump-version.sh','distribute.sh'], 'Docs/build':['*']},
    install_requires = ['pyfits>=2.4','numpy>=1.6','scipy>=0.9','matplotlib>=1.1','AstroObject>=0.2.5'],
    dependency_links = ['https://github.com/alexrudy/AstroObject/tags'],
    test_suite = 'Tests',
    entry_points = {
        'console_scripts' : ['SEDMsim = SED.Simulator:run',],
    },
)
