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
    scripts=['bin/mgc3.py','bin/mgc3_get_pole_stars.py','bin/mgc3_plot_polemaps.py','bin/peakdetect_mgc3_polemaps.py'],
    url='https://github.com/cmateu/PyMGC3',
    license='LICENSE.txt',
    description='PyMGC3 codes',
    long_description=open('README.rst').read(),
    install_requires=[
      "numpy",
      "scipy"
    ],
)
