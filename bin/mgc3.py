#!/usr/bin/env python
from numpy import array
import sys
import scipy
import mgc3_lib as mgc3_lib
import argparse
#############################################################################
#Copyright (c) 2013 - 2014, Cecilia Mateu
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without 
#modification, are permitted provided that the following conditions are met:
#
#   Redistributions of source code must retain the above copyright notice, 
#      this list of conditions and the following disclaimer.
#   Redistributions in binary form must reproduce the above copyright notice, 
#      this list of conditions and the following disclaimer in the 
#      documentation and/or other materials provided with the distribution.
#   The name of the author may not be used to endorse or promote products 
#      derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
#BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
#OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
#AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
#WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#POSSIBILITY OF SUCH DAMAGE.
#############################################################################
#=======MAIN=============

'''
NAME:
   mgc3
PURPOSE:
   Produce pole count map using the mGC3 method from Mateu et al. 2011 (MNRAS, 415, 214-2249
INPUT:
   parameter_file - Parameter file
   data_file      - Observations file
   output_ext     - Prefix to be used for output filename before standard extension ".mgc3.cts"
OUTPUT:
   phi_pole,theta_pole,n_mgc3  
   This data are printed on output file data_file.output_ext
HISTORY:
   2013-11-21 - Uses mgc3_lib (v4) (nmgc3 pole counts added)
		Input gzip catalogue supported
   2013-09-20 - Added gc3 pole counts attribute to pole-grid object
   2013-09-xx - Written. Cecilia Mateu (CIDA,IA-UNAM Ensenada)
'''

#Parse command line arguments

__version__ = '4.0'
__docformat__ = "reredtext en"
__what__= sys.argv[0]+": This program computes mGC3 pole counts for the input catalogue"
#
parser = argparse.ArgumentParser(description=__what__)
#

parser.add_argument('parfile',metavar='parameter_file',help='Input catalogue parameter file',action='store',nargs=1)
parser.add_argument('infile',metavar='data_file',help='Input catalogue file',nargs=1,action='store')
parser.add_argument('ext_prefix',metavar='outfile_ext',help='Output file prefix (output will be infile.outfile_ext.mgc3.cts)',nargs=1,action='store')
parser.add_argument('-farea',help='For non-allsky surveys, compute area fraction correction', action='store_true')
#parser.add_argument('-ppar','--print_parf',help='Print sample parameter file mgc3.par and exit', action='store_true')
parser.add_argument('-ppar','--print_parf',help='Print sample parameter file mgc3.par and exit',action=mgc3_lib.print_parfile_action,nargs=0)
parser.add_argument('-v',"--version", help="Print program version", action="version",version='Version '+__version__)
parser.add_argument('-doc',"--documentation", help="Print short program description", action="version",version=__what__)

args = parser.parse_args()
parfile, filename, ext_prefix = args.parfile[0], args.infile[0], args.ext_prefix[0]

#If asked for by user, print sample parameter file
#if args.print_parf: 
#  mgc3_lib.print_sample_parfile()
#  sys.exit('Exiting.')

#Parse parameter file
print 'Reading Parameter file %s ...' % (parfile)
survey_pars=mgc3_lib.parse_pars(parfile)

print 'Reading input file %s ...' % (filename)
obsdata,filename=mgc3_lib.read_inputcat_for_mgc3(filename,pardic=survey_pars)
print 'Input file shape (rows,cols): ', obsdata.shape

print 'Initializing pole grid...'
mygrid=mgc3_lib.pole_grid(poles=survey_pars['grid_step'],pole_grid_dic=survey_pars)

print 'Computing pole counts...'
mygrid.mgc3(obsdata,pars=survey_pars)

if args.farea:
  print 'Initializing auxiliary pole grid for farea computation...'
  mygrid_foot=mgc3_lib.pole_grid(poles=survey_pars['grid_step'])
  foot_survey,foot_survey_pars=mygrid_foot.get_uniform_survey_footprint(obsdata,pars=survey_pars)
  mygrid_foot.mgc3(foot_survey,pars=foot_survey_pars)
else:
  mygrid_foot=mygrid
  mygrid_foot.farea=mygrid_foot.farea*0.+1.  

outfilename=filename.replace('.dat','')+'.'+ext_prefix+'.mgc3.cts'
print 'Printing output file %s ...' % (outfilename)
ofile=open(outfilename,'w')
ofile.write("#%9s %10s %10s %10s %10s %10s %10s\n" % ("phi","theta","np_mgc3gal","np_mgc3hel",'np_gc3gal',"np_ngc3gal",'farea'))
scipy.savetxt(ofile,array([mygrid.l,mygrid.b,mygrid.np_mgc3,mygrid.mgc3hel,mygrid.np_gc3,mygrid.np_ngc3,mygrid_foot.farea]).T,fmt='%10.3f %10.3f %10d %10d %10d %10d %10.4f')
print 'Done'
