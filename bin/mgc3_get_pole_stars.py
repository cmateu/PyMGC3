#!/usr/bin/env python
import sys
import scipy
import mgc3_lib as mgc3_lib
import myutils
import numpy as np
import argparse

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
#parser.add_argument('ext_prefix',metavar='outfile_ext',help='Output file prefix (output will be infile.outfile_ext.mgc3.cts)',nargs=1,action='store')
parser.add_argument('-ext',metavar='outfile_ext',help='Output suffix [optional]. If given output will be infile.outfile_ext.mgc3.pst',action='store',default=['',],nargs=1)
parser.add_argument('-ppar','--print_parf',help='Print sample parameter file mgc3.par and exit', action='store_true')
parser.add_argument('-m',help='Select stars using mGC3/nGC3/GC3 method criteria. Default is mGC3', action='store',default='mGC3',choices=['mGC3','nGC3','GC3'])
parser.add_argument('-v',"--version", help="Print program version", action="version",version='Version '+__version__)
parser.add_argument('-doc',"--documentation", help="Print short program description", action="version",version=__what__)

args = parser.parse_args()

parfile, filename, polelistname, ext_prefix = args.parfile[0], args.infile[0], args.polelist[0], args.ext[0]

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
ext='%s.%s.pst' % (ext_prefix,args.m.lower())
outfilename=filename.replace('.dat','')+ext
print 'Printing output file %s ...' % (outfilename)
#Open output file and print the exact same header the input file had, with an extra column
outfile=open(outfilename,'w')
outfile.write('#Stars selected according to %s criteria\n' % (args.m))
scipy.savetxt(outfile,head,fmt='%s')
#The masks associated to each pole will be combined with OR, this way, each star can only be printed once,
#even if associated to more than one pole
indep_mask=np.zeros(obsdata[:,0].size,dtype=bool)
indep_pole_ids=-1*np.ones(obsdata[:,0].size)
for id_pole,phi_pole,theta_pole in polelist:
  mygrid=mgc3_lib.pole_grid(poles=[phi_pole,theta_pole])
  cat_mask=mygrid.mgc3_allobs_one_pole(obsdata,pars=survey_pars,return_mask=args.m)
  print '   stars associated to pole %s: %d' % (id_pole,np.sum(1*cat_mask))
  #combine mask with OR
  indep_mask=indep_mask | cat_mask
  #label stars added to the mask in this step with current poleid
  indep_pole_ids[cat_mask]=id_pole
  #pole_ids=id_pole*np.ones_like(obsdata[cat_mask,0])

#Printing is now done after finishing the loop
print_data=obsdata[indep_mask,:].T
indep_pole_ids=indep_pole_ids[indep_mask]
scipy.savetxt(outfile,np.vstack([print_data,indep_pole_ids]).T)
