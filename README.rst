PyMGC3 
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
uses positional information to search for great-circle-cell
structures; mGC3 makes use of full 6D data and 
nGC3 uses positional and proper motion data.

----------

**LATEST FEATURES:**

2019/07 - Solar position and distance now printed in cts file header. Default values updated.

2019/05 - New option -arep added to mgc3_get_pole_stars.py - Allows repetitions when printing out stars associated to a set of poles (a star can appear multiple times) 

2019/04 - Glad to announce: *PyMGC3 is now compatible with Python 3* !!!

2017/03 - New utility code combine_mgc3_maps.py combines two pcms using pole counts from different methods. Useful, e.g., to combine GC3 + nGC3 when only some stars have proper motion data.

2017/02 - Standard deviation of background counts now computed locally in annulus around each pixel. New option -npixann added to set annulus radius. Option npixmin added to control minimum number of pixels required for valid detections.

2016/02 - Heliocentric GC3 counts added (np_gc3hel column)

2016/02 - Input list option added to mgc3.py (use -l flag to take input as a list of input files to run mgc3.py). Normalization option added to join_mgc3_maps.py (normalizes each input cts file's max counts and then combines them).   

2015/03 - Unsharp mask option (and new cmd-line parameters) added to peak detection 
code (peakdetect_mgc3_polemaps.py). Peak detection can now be done in unsharp-masked PCM.

2015/03 - New python scripts added for plotting (mgc3_plot_polemaps.py) and automatically detecting 
peaks (peakdetect_mgc3_polemaps.py) in pole count maps.

----------


**REQUIREMENTS**

- Python modules required are NUMPY and SCIPY. MATPLOTLIB and BASEMAP are needed for plotting utilities.
- This programs makes use of the coordinate transformation library
  bovy_coords.py from the `galpy <https://github.com/jobovy/galpy>`__ 
  package by Jo Bovy (2015, in prep.). It is supplied with this bundle.
- The peak detection utility peakdetect_mgc3_polemaps.py uses the
  Starlink implementation of the Fellwalker code by `Berry 2014 <http://arxiv.org/abs/1411.6267v1>`__,
  assuming it is installed at /star-2014A. The code is *not* supplied
  with this bundle, but its publicly available at the `Starlink website <http://starlink.jach.hawaii.edu>`__.

**FILES PROVIDED**

- Executable programs
   * mgc3.py
   * mgc3_get_pole_stars.py
   * mgc3_plot_polemaps.py
   * peakdetect_mgc3_polemaps.py
   * mgc3_peakstars_plot.py  
   * join_mgc3_maps.py       
   * combine_mgc3_maps.py
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

In a terminal, run the following command::

    sudo python setup.py install

Source your .cshrc

If you do not have root access, you can install in the custom directory path_to_dir.
First, add the directory's path path_to_dir and path_to_dir/lib/python2.7/site-packages/ 
to the PYTHONPATH variable in your .cshrc/.bashrc file and source it. Then install using the --prefix option::

    python setup.py install --prefix=path_to_dir

Add path_to_dir/bin to your PATH in your .csrhc or .bashrc file.

Quick Guide
-----------

If you just want to start running the programs right away, follow these quick-and-dirty recipe
using the example data provided. For more details on each of the Python scripts,
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

    mgc3_plot_polemaps.py example_data.test01.mgc3.cts


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

    usage: mgc3.py [-h] [-l] [-farea] [-ppar] [-v] [-doc]
                   parameter_file data_file outfile_ext
    mgc3.py: error: too few arguments

Run with -h or --help argument for full help, like this::

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

Name of the input catalogue file. Assumed to be ascii format (.gzip supported), with comments preceeded by #.
This file can also be read as a list of input catalogue filenames with the -l option.

*ext_prefix*

The output file returned by mgc3.py will be called data_file.ext_prefix.mgc3.cts. 
It will contain (phi,theta) and pole counts np_mgc3_gal (MGC3), 
np_gc3gal (GC3), np_ngc3gal (nGC3) for an uniform pole grid with a step 
given by grid_step. It also contains mgc3 heliocentric (np_mgc3hel) pole counts, these
are useful for experimentation sometimes. Note the pole grid covers one hemisphere, 
as the information from the other hemisphere is redundant.

Note: a program for plotting maxima in pole count maps (mgc3_plot_polemaps.py)
is provided with the PyMGC3 bundle. You can also quickly plot using Topcat 
(`<http://www.star.bris.ac.uk/~mbt/topcat/>`_), Aitoff or Sin projections are recommended.

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
it uses mGC3 criteria, but any of the mgc3-family methods (mGC3/nGC3/GC3/mGC3hel/GC3hel) can be used to 
select stars associated to each of the poles in the list.

Run without arguments for a short help message to explain inputs and optional arguments::

    mgc3_get_pole_stars.py
    usage: mgc3_get_pole_stars.py [-h] [-ext outfile_ext] [-ppar]
                                  [-m {mGC3,nGC3,GC3,mGC3hel,GC3hel}] [-v] [-doc]
                                  parameter_file data_file pole_list
    mgc3_get_pole_stars.py: error: too few arguments


Run with -h or --help for full help::

    mgc3_get_pole_stars.py -h

**OUTPUT:**

The output file infile.mgc3.pst is identical to the input catalogue, but including only stars associated 
with the given poles plus an additional column at the end indicating the pole_ID for the pole
each star is associated with.

**EXAMPLE:**

Use the pole count map and pole list examples as the input for this program::

    mgc3_get_pole_stars.py example_data.par example_data.dat example_data.mgc3.pls my_test

The output file will be example_data.my_test.mgc3.dat. Try running with the -m nGC3 and -m GC3 flags
to get \*.ngc3.dat and \*.gc3.dat outputs.

The extension of the input .pls file must be consistent with the method selected to extract the stars,
the program will exit with a warning if this is not the case to be safe.

Program mgc3_plot_polemaps.py
-----------------------------

**DESCRIPTION:**

This program plots mGC3-family (mGC3,nGC3,GC3,mGC3hel,GC3hel) pole count maps in different projections. 

**SYNTAX:**

The required command line arguments are:

*polecounts_file*

Pole count map file (.cts). Assumes the default output format from the mgc3.py program.
If the -l flag is set, this is assumed to be a list of pole count map files

*Optional arguments*

The program has several optional keywords and flags to customize the output plot, for a full list
and details run with the -h or --help flag::

        usage: mgc3_plot_polemaps.py [-h] [-l]
                                     [-m {mGC3,nGC3,GC3,GC3hel,mGC3hel,smooth,usharpc,usharpn}]
                                     [-f {png,eps,pdf}] [-proj {npaeqd,ortho,moll}]
                                     [-log] [-lon0 LON0] [-lat0 LAT0] [-dlat DLAT]
                                     [-dlon DLON] [-latmax LATMAX] [-mlab] [-mlabr]
                                     [-plab] [-vmin VMIN] [-vmax VMAX] [-ms MS] [-c]
                                     [-t] [-s] [-title TITLE] [-pls PLSFILE]
                                     [-al ALPHA] [-ff FFONTS] [-flab FLABELS]
                                     [-fcirc FCIRC]
                                     [-cmap {sron,gray,gray_r,viridis,inferno}]
                                     [-ext outfile_ext]
                                     infile
        
        Plot mGC3/nGC3/GC3 pole count maps
        
        positional arguments:
          infile                Input file containing pole count maps (\*.cts file)
        
        optional arguments:
          -h, --help            show this help message and exit
          -l, --llist           Take infile as list of mgc3.cts files
          -m {mGC3,nGC3,GC3,GC3hel,mGC3hel,smooth,usharpc,usharpn}
                                Plot pole count map. Default is mGC3
          -f {png,eps,pdf}, --fig {png,eps,pdf}
                                Output plot type png/eps. Default is png
          -proj {npaeqd,ortho,moll}
                                Projection npaeqd/ortho/mollweide. Default is npaeqd
          -log                  Plot pole-count map in log-scale
          -lon0 LON0            Longitude for Y-axis. Default is 0.
          -lat0 LAT0            Bounding latitude for plot. Default is 90.
          -dlat DLAT            Spacing between parallels. Default is 20.
          -dlon DLON            Spacing between meridians. Default is 20.
          -latmax LATMAX        Max latitude upto which meridians are drawn. Default
                                is 80.
          -mlab, --merlabels    Show meridian labels. Default is False
          -mlabr, --merlabelsr  Show meridian labels (right axes). Default is False
          -plab, --parlabels    Show parallel labels. Default is False
          -vmin VMIN            Min counts for color-scale. Default is min(cts)
          -vmax VMAX            Max counts for color-scale. Default is max(cts)
          -ms MS                Marker size. Default: 15/40 for npaeqd/ortho.
          -c, --contour         Plot pole-count contour map instead of raw grid.
          -t, --twohemispheres  Plot both hemispheres in pole-count map.
          -s, --show            Show plot in window. Default is False
          -title TITLE          Plot title
          -pls PLSFILE          Overplot poles from peakdetect output file (.pls)
          -al ALPHA, --alpha ALPHA
                                Clump transparency. Default 0.4
          -ff FFONTS, --ffonts FFONTS
                                Increase size tick and axes labels by factor ff.
                                Default 1.
          -flab FLABELS, --flabels FLABELS
                                Increase size of peak labels by factor flab. Default
                                1.
          -fcirc FCIRC, --fcirc FCIRC
                                Increase size of peak markers by factor fcirc. Default
                                1.
          -cmap {sron,gray,gray_r,viridis,inferno}
                                Choose color map. Default is sron
          -ext outfile_ext      Output suffix [optional]. If given output will be
                                infile.outfile_ext.mgc3.pst       

**EXAMPLES:**

Use the example data to produce a pole counts file with::

  mgc3.py example_data.par example_data.dat test02

The following example plots the resulting map for the nGC3 pole counts, using the Mollweide projection, with meridians every 30 deg and paralles every 20deg. The -t flag forces both hemispheres to be plotted in the map. The output is saved in pdf format:: 

  mgc3_plot_polemaps.py example_data.test02.mgc3.cts -m nGC3 -dlat 30 -dlon 20 -proj moll -t -f pdf
 
The output figure is called example_data.test02.mgc3.moll.r.pdf.  

Selection the ortho projection produces a figure with the map as seen from lon0 and lon0+180deg to ensure the whole map is visible::

  mgc3_plot_polemaps.py example_data.test02.mgc3.cts -m GC3 -f pdf -dlat 30 -dlon 20 
                         -proj ortho -lon0 65

The output figure is called example_data.test02.mgc3.ortho.r.pdf. 

Pole count contour plots can be plotted with the -c option::

  mgc3_plot_polemaps.py example_data.test02.mgc3.cts -m nGC3 -f png -dlat 30 -dlon 20 -c

The output figure is called example_data.test02.mgc3.npa.c.png. Note: the -c option works
only in the npaeqd projection for now.


Program peakdetect_mgc3_polemaps.py
-----------------------------------

**DESCRIPTION:**

This program detects peaks in pole-count maps after unsharp masking and plots the pole count map
indicating the peaks found. It needs the Fellwalker code to run (Berry 2014).

**SYNTAX:**

The only required argument is the pole-count file (or list when using the -l option). 

Run with -h for a full list of options::


  peakdetect_mgc3_polemaps.py -h

Run with -nc for plotting only:: 

  peakdetect_mgc3_polemaps.py example_data.test02.mgc3.cts -nc 

Most plotting options available are the same as for mgc3_plot_polemaps.py. Two 
ways are available to select the minimum peak height threshold value::

  peakdetect_mgc3_polemaps.py example_data.test02.mgc3.cts -frms 5

The option -frms 5 means the peaks must have a height >5*RMS, where RMS is
the root mean squared deviation of the pole counts. This threshold can
also be defined as a fraction of the maximum counts in the map with 
the -ffrac option:: 

  peakdetect_mgc3_polemaps.py example_data.test02.mgc3.cts -ffrac 0.6

In this case, peaks must be at least 0.6*max_counts to be saved. 

Program join_mgc3_maps.py
-----------------------------------

**DESCRIPTION:**

This utility program sums pole counts in a list of pole-count maps.

**SYNTAX:**

The required arguments are a list of pole-count map files (.mgc3.cts) and a name
for the output file::

  join_mgc3_maps.py  infile_list outfilename

Its highly recommended to use the .mgc3.cts extension for the output file, for 
consistency with the rest of PyMGC3 programs. 

The -n option normalizes the different pole counts (GC3,nGC3,mGC3, etc.) 
in each of the input pole-count maps before adding them up.

**INPUTS AND OUTPUTS:**

Input files are assumed to have the same format as mgc3.py outputs. The output
file will have the same format as well.

When the -n option is used, an extra set of output files (e.g. outfilename.mgc3.pls)
is produced listing the coordinates and counts for the maximum used for normalization
of each input file.

Program combine_mgc3_maps.py
-----------------------------------

**DESCRIPTION:**

This utility program combines pole counts from different methods for two input pole count maps.

**SYNTAX:**

The required arguments are the input file names, method counts to be added and output file name::

  combine_mgc3_maps.py infile1 method1 infile2 method2

Choices for method1 and method2 are: {mGC3,nGC3,GC3,mGC3hel,GC3hel}

Its highly recommended to use the .mgc3.cts extension for the output file, for 
consistency with the rest of PyMGC3 programs. 

**INPUTS AND OUTPUTS:**

Input files are assumed to have the same format as mgc3.py outputs. The output
file will have the same format as well.

Attribution
-----------

Cecilia Mateu - cmateu at cida.gob.ve

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
