#!/usr/bin/env python
import sys
import scipy
import mgc3_lib as mgc3_lib
import myutils
import numpy as np
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

#Parse command line arguments
__version__ = '3.0'
__docformat__ = "reredtext en"
__what__= sys.argv[0]+": This program extracts stars associated to a list of mgc3 poles. Author: C. Mateu - 2013"
#
parser = argparse.ArgumentParser(description=__what__)
#

parser.add_argument('parfile',metavar='parameter_file',help='Input catalogue parameter file',action='store',nargs=1)
parser.add_argument('infile',metavar='data_file',help='Input catalogue file',nargs=1,action='store')
parser.add_argument('polelist',metavar='pole_list',help='List of pole coordinates [ID phi_pole theta_pole]',nargs=1,action='store')
parser.add_argument('ext_prefix',metavar='outfile_ext',help='Output file prefix (output will be infile.outfile_ext.mgc3.cts)',nargs=1,action='store')
parser.add_argument('-ppar','--print_parf',help='Print sample parameter file mgc3.par and exit', action='store_true')
parser.add_argument('-m',help='Select stars using mGC3/nGC3/GC3 method criteria. Default is mGC3', action='store',default='mGC3',choices=['mGC3','nGC3','GC3'])
parser.add_argument('-v',"--version", help="Print program version", action="version",version='Version '+__version__)
parser.add_argument('-doc',"--documentation", help="Print short program description", action="version",version=__what__)

args = parser.parse_args()

parfile, filename, polelistname, ext_prefix = args.parfile[0], args.infile[0], args.polelist[0], args.ext_prefix[0]

#If asked for by user, print sample parameter file
if args.print_parf:
  mgc3_lib.print_sample_parfile()
  sys.exit('Exiting.')

#Parse parameter file
print 'Reading Parameter file %s ...' % (parfile)
survey_pars=mgc3_lib.parse_pars(parfile)

print 'Reading input file %s ...' % (filename)
obsdata = scipy.genfromtxt(filename,comments='#')
print 'Input file shape (rows,cols): ', obsdata.shape

#Save input file header
head = myutils.get_header_line(filename)
head[-1]=head[-1]+' IDpole' 

print 'Reading pole list file %s ... ' % (polelistname)
polelist=scipy.genfromtxt(polelistname,comments='#',usecols=(0,1,2))
if np.ndim(polelist)==1: polelist=np.array([polelist,])

#Cicle over pole list, one by one
ext='.%s.%s.pst' % (ext_prefix,args.m.lower())
outfilename=filename.replace('.dat','')+ext
print 'Printing output file %s ...' % (outfilename)
#Open output file and print the exact same header the input file had, with an extra column
outfile=open(outfilename,'w')
outfile.write('#Stars selected according to %s criteria\n' % (args.m))
scipy.savetxt(outfile,head,fmt='%s')
for id_pole,phi_pole,theta_pole in polelist:
  mygrid=mgc3_lib.pole_grid(poles=[phi_pole,theta_pole])
  cat_mask=mygrid.mgc3_allobs_one_pole(obsdata,pars=survey_pars,return_mask=args.m)
  print '   stars associated to pole: ',sum(1*cat_mask)
  pole_ids=id_pole*np.ones_like(obsdata[cat_mask,0])
  print_data=obsdata[cat_mask,:].T
  scipy.savetxt(outfile,np.vstack([print_data,pole_ids]).T)
