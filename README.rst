PyMGC3 Package
======

**DESCRIPTION:**

PyMGC3 is a Python toolkit to apply the Modified Great Circle 
Cell Counts (mGC3) method, devised by Mateu et al. (2011) to 
search for tidal streams in the Galactic Halo. 

The current distribution computes pole count maps using 
the full mGC3/nGC3/GC3 family of methods described 
in Mateu et al. 2011 (MNRAS, 415, 214-224) and 
Abedi et al. 2014 (MNRAS, submitted). Briefly, 
the original GC3 method developed by Johnston et al. (1996) 
uses positional information to search for 'great-circle-cell
structures'; mGC3 makes use of full 6D data and 
nGC3 uses positional and proper motion data.


**REQUIREMENTS**

- Python modules required are NUMPY and SCIPY.
- This programs makes use of the coordinate transformation library
  bovy_coords.py by Jo Bovy (2011). It is supplied with this bundle.

**FILES PROVIDED**

- Executable programs
   * mgc3.py
   * mgc3_get_pole_stars.py
- Documentation
   * README.pst
- Libraries
   * bovy_coords.py
   * mgc3_lib.py
- Example data
   * example_data.dat
   * example_data.par
   * example_data.pls
   * example_cts_output.png
   * example_pst_output.png

**INSTALLATION**

In a terminal run the following command::

    sudo python setup.py install

If you do not have root access, install in a custom directory using the --prefix option::

    python setup.py install --prefix=path_to_dir

After installing, add path_to_dir/PyMGC3/bin to your PATH in your .csrhc or .bashrc file.
Also add path_to_dir/PyMGC3/bin and path_to_dir/PyMGC3/ to PYTHONPATH also in your .cshrc/.bashrc file.

Quick Guide
-----------

If you just want to start running the programs right away, follow these quick-and-dirty recipe
using the example data provided. For more details on the Python scripts,
goto the next section in this README.

Go to the examples directory::

    cd examples

Print sample parameter file::

    mgc3.py -ppar

Run mgc3 to get pole count maps::

    mgc3.py mgc3_sample.par example_data.dat test01

Use the pole list provided with the example to extract stars associated
to the maxima in the pole count map::

    mgc3_get_pole_stars.py mgc3_sample.par example_data.dat example_data.mgc3.pls test01
   

Program mgc3.py
---------------


**DESCRIPTION:**

This program is used to apply the mGC3 method to an observations catalogue
supplied by the user. The catalogue must have all 6D information for each
star, i.e. position angles and proper motions (galactic or equatorial),
parallax and radial velocity. The program returns mGC3, nGC3 and GC3 pole
counts. For details on the mGC3/nGC3/GC3 methods see Mateu et al. 2011
and Abedi et al. 2014.

**SYNTAX:**

The required command line arguments are:

*parameter_file*: the name of the parameter file to be used

*data_file*: the name of catalogue/data file to be used

*outfile_extension*: an extension to be used for the output file

Running the mgc3.py without any arguments will provide a short description
of the required syntax and ask the user whether a sample parameter file
should be printed::

    mgc3.py

    usage: mgc3.py [-h] [-farea] [-ppar] [-v] [-doc]
               parameter_file data_file outfile_extension
    mgc3.py: error: too few arguments

Run with -h or --help argument for full help like this::

    mgc3.py -h

Run with -ppar flag to print a sample parameter file::

    mgc3.py -ppar

the output file will be mgc3_sample.par

**INPUTS AND OUTPUTS:**

*parameter_file*

The parameter file indicates the structure of the input catalogue,
as well as the values to be used for mGC3 parameters. Each parameter
is explained briefly by a comment in the sample parameter file header. 

*data_file*

Name of the input catalogue file. Assumed to be ascii format, with comments preceeded by #.

*ext_prefix*

The output file returned by mgc3.py will be called data_file.ext_prefix.mgc3.cts. 
It will contain (phi,theta) and pole counts np_mgc3_gal (MGC3), 
np_gc3gal (GC3), np_ngc3gal (nGC3) for an uniform pole grid with a step 
given by grid_step. It also contains mgc3 heliocentric (np_mgc3hel) pole counts, these
are useful for experimentation sometimes. Note the pole grid covers one hemisphere, 
as the information from the other hemisphere is redundant.

Note: a program for plotting and detecting maxima in pole count maps will be provided
with the mgc3 bundle in an upcoming version. In the mean time you can quickly
plot using Topcat (`<http://www.star.bris.ac.uk/~mbt/topcat/>`_), it is recommended to
use an Aitoff or Sin projection.

**EXAMPLE:**

To test mgc3.py, use the provided example data, running mgc3.py with this command line::

    mgc3.py example_data.par example_data.dat  my_test

The output file will be example_data.my_test.mgc3.cts. The output pole maps 
should look like those shown on example_output.png when plotted with Topcat
using either 3D spherical mode (left) or aitoff map mode (right). See below
on how to plot the maps and detect maxima.

Program mgc3_get_pole_stars.py
------------------------------

**DESCRIPTION:**

This program extracts stars associated to poles given in an input list. By default
it uses mGC3 criteria, but any of the three methods (mGC3/nGC3/GC3) can be used to 
select stars associated to each of the poles in the list.

Run without arguments for a short help message to explain inputs and optional arguments::

    get_mgc3pole_stars.py
    usage: mgc3_get_pole_stars.py [-h] [-ppar] [-m {mGC3,nGC3,GC3}] [-v] [-doc]
                                  parameter_file data_file outfile_ext pole_list

    mgc3_get_pole_stars.py: error: too few arguments

Run with -h or --help for full help::

    get_mgc3pole_stars.py -h

**OUTPUT:**

The output file infile.mgc3.pst is identical to the input catalogue, but including only stars associated 
with the given poles and with an additional column at the end indicating the pole_ID for the pole
each star is associated with.

**EXAMPLE:**

Use the pole count map and pole list examples as the input for this program::

    mgc3_get_pole_stars.py example_data.par example_data.dat example_data.mgc3.pls my_test

The output file will be example_data.my_test.mgc3.dat. Try running with the -m nGC3 and -m GC3 flags
to get \*.ngc3.dat and \*.gc3.dat outputs.

Attribution
-----------

Cecilia Mateu - cmateu at astrosen.unam.mx

If you have used this code in your research, please let me know and consider acknowledging this package.

License
-------

Copyright (c) 2013-2014 Cecilia Mateu

PyMGC3 is open source and free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see `<http://www.gnu.org/licenses/>`_.
