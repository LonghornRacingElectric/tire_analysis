#!/usr/bin/env python

from setuptools import find_packages, setup

setup(name='LHR_tire_toolkit',
      version='1.0',
      description='Toolkit fot MF5.2 Tire Models',
      author='Robert Horvath',
      author_email='rhorvath@utexas.edu',
      url='https://github.com/LonghornRacingElectric/tire_analysis',
      packages=['tire_model', 'tire_model/file_processing', 'tire_model/MF52_calculations'],
      py_modules=['tire_model/MF52']
      # packages=find_packages(exclude=("assets", "__pycache__",)),
     )