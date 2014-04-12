import os
import sys
import re

try:
  from setuptools import setup
except ImportError:
  from distutils.core import setup


setup(
    name='PyMGC3',
    version='0.1.0',
    author='C. Mateu',
    author_email='cmateu@astrosen.unam.mx',
    packages=['bovy_coords','mgc3_lib','myutils'],
    scripts=['bin/mgc3.py','bin/mgc3_get_pole_stars.py'],
    url='http://www.cida.ve/~cmateu/Gaia',
    license='LICENSE.txt',
    description='MGC3 codes',
    long_description=open('README.rst').read(),
    install_requires=[
      "numpy",
      "scipy"
    ],
)
