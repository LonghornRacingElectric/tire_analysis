#!/usr/bin/env python

from setuptools import find_packages, setup

setup(name='LHR_tire_toolkit',
      version='1.0',
      description='Toolkit fot MF5.2 Tire Models',
      author='Robert Horvath',
      author_email='rhorvath@utexas.edu',
      url='https://github.com/LonghornRacingElectric/tire_analysis/tree/LHR_tire_toolkit_package',
      packages=['LHR_tire_toolkit'],
      install_requires=[
          'numpy',
          'scipy',
          'pandas'
      ]
     )