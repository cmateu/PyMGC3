PyMGC3 Package
======

**DESCRIPTION:**

PyMGC3 is a Python toolkit to apply the Modified Great Circle 
Cell Counts (mGC3) method, devised by `Mateu et al. (2011) <http://adsabs.harvard.edu/abs/2011MNRAS.415..214M>`__ 
to search for tidal streams in the Galactic Halo. 

The current distribution computes pole count maps using 
the full mGC3/nGC3/GC3 family of methods described 
in `Mateu et al. (2011) <http://adsabs.harvard.edu/abs/2011MNRAS.415..214M>`__ and 
`Abedi et al. (2014) <http://adsabs.harvard.edu/abs/2014MNRAS.442.3627A>`__. Briefly, 
the original GC3 method developed by `Johnston et al. (1996) <http://adsabs.harvard.edu/abs/1996ApJ...465..278J>`__
uses positional information to search for 'great-circle-cell
structures'; mGC3 makes use of full 6D data and 
nGC3 uses positional and proper motion data.


**REQUIREMENTS**

- Python modules required are NUMPY and SCIPY.
- This programs makes use of the coordinate transformation library
  bovy_coords.py from the `galpy <https://github.com/jobovy/galpy>`__ 
  package by Jo Bovy (2015, in prep.). It is supplied with this bundle.

**FILES PROVIDED**

- Executable programs
   * mgc3.py
   * mgc3_get_pole_stars.py
   * plot_mgc3_polemaps.py
- Documentation
   * README.rst
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
go to the next section in this README.

Go to the examples directory::

    cd examples

Print sample parameter file::

    mgc3.py -ppar

Run mgc3 to get pole count maps::

    mgc3.py mgc3_sample.par example_data.dat test01

Use the pole list provided with the example to extract stars associated
to the maxima in the pole count map::

    mgc3_get_pole_stars.py mgc3_sample.par example_data.dat example_data.mgc3.pls test01
   
Plot the mGC3 pole count map in a north-polar-azimuthal equidistant projection::

    plot_mgc3_polemaps.py example_data.mgc3.cts


Program mgc3.py
---------------


**DESCRIPTION:**

This program is used to apply the mGC3 method to an observations catalogue
supplied by the user. The catalogue must have all 6D information for each
star, i.e. position angles and proper motions (galactic or equatorial),
parallax and radial velocity. The program returns mGC3, nGC3 and GC3 pole
counts. For a detailed explanation on the mGC3/nGC3/GC3 methods see 
`Mateu et al. (2011) <http://adsabs.harvard.edu/abs/2011MNRAS.415..214M>`__ and
`Abedi et al. (2014) <http://adsabs.harvard.edu/abs/2014MNRAS.442.3627A>`__.

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

Program plot_mgc3_polemaps.py
-----------------------------

**DESCRIPTION:**

This program plots mGC3/nGC3/GC3 pole count maps in different projections. 

**SYNTAX:**

The required command line arguments are:

*polecounts_file*

Pole count map file (.cts). Assumes the default output format from the mgc3.py program

*Optional arguments*

The program has several optional keywords and flags to customize the output plot, for a full list
and details run with the -h or --help flag::

    plot_mgc3_polemaps.py -h

    usage: plot_mgc3_polemaps.py [-h] [-l] [-m {mGC3,nGC3,GC3}] [-f {png,eps,pdf}]
                                 [-proj {npaeqd,ortho,mollweide}] [-lon0 LON0]
                                 [-lat0 LAT0] [-dlat DLAT] [-dlon DLON] [-ms MS]
                                 [-c] [-t] [-s]
                                 infile
    
    positional arguments:
      infile                Input file containing pole count maps (*.cts file)
    
    optional arguments:
      -h, --help            show this help message and exit
      -l, --llist           Take infile as list of mgc3.cts files
      -m {mGC3,nGC3,GC3}    Plot mGC3/nGC3/GC3 pole count map. Default is mGC3
      -f {png,eps,pdf}, --fig {png,eps,pdf}
                            Output plot type png/eps. Default is png
      -proj {npaeqd,ortho,moll}
                            Projection npaeqd/ortho/mollweide. Default is npaeqd
      -lon0 LON0            Longitude for Y-axis. Default is 0.
      -lat0 LAT0            Bounding latitude for plot. Default is 90.
      -dlat DLAT            Spacing between parallels. Default is 30.
      -dlon DLON            Spacing between meridians. Default is 30.
      -ms MS                Marker size. Default: 90/40 for npaeqd/ortho.
      -c, --contour         Plot pole-count contour map instead of raw grid.
      -t, --twohemispheres  Plot both hemispheres in pole-count map.
      -s, --show            Show plot in window. Default is False


**EXAMPLES:**

Use the example data to produce a pole counts file with::

  mgc3.py example_data.par example_data.dat test02

The following example plots the resulting map for the nGC3 pole counts, using the Mollweide projection, with meridians every 30 deg and paralles every 20deg. The -t flag forces both hemispheres to be plotted in the map. The output is saved in pdf format:: 

  plot_mgc3_polemaps.py example_data.test02.mgc3.cts -m nGC3 -dlat 30 -dlon 20 -proj moll -t -f pdf
 
The output figure is called example_data.test02.mgc3.moll.r.pdf.  

Selection the ortho projection produces a figure with the map as seen from lon0 and lon0+180deg to ensure the whole map is visible::

  plot_mgc3_polemaps.py example_data.test02.mgc3.cts -m GC3 -f pdf -dlat 30 -dlon 20 
                         -proj ortho -lon0 65

The output figure is called example_data.test02.mgc3.ortho.r.pdf. 

Pole count contour plots can be plotted with the -c option::

  plot_mgc3_polemaps.py example_data.test02.mgc3.cts -m nGC3 -f png -dlat 30 -dlon 20 -c

The output figure is called example_data.test02.mgc3.npa.c.png. Note: the -c option is working 
only in the npaeqd projection for now.

Attribution
-----------

Cecilia Mateu - cmateu at astrosen.unam.mx

If you have used this code in your research, please let me know and consider acknowledging this package.

License
-------

Copyright (c) 2013-2014 Cecilia Mateu

PyMGC3 is open source and free software: 
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. The name of the author may not be used to endorse or promote
products derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
