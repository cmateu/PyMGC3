#!/usr/bin/env python
import sys
import scipy
import mgc3_lib as mgc3_lib
import myutils
import numpy as np
import argparse
import gzip 

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
parser.add_argument('-ext',metavar='outfile_ext',help='Output suffix [optional]. If given output will be infile.outfile_ext.mgc3.pst',action='store',default=['',],nargs=1)
parser.add_argument('-ppar','--print_parf',help='Print sample parameter file mgc3.par and exit', action='store_true')
parser.add_argument('-m',help='Select stars using mGC3/nGC3/GC3/mGC3hel/GC3hel method criteria. Default is mGC3', action='store',default='mGC3',choices=['mGC3','nGC3','GC3','mGC3hel','GC3hel'])
parser.add_argument('-arep','--allow_repeats',help='Force printing each star for each pole, allowing a star to be assigned to multiple poles (default behaviour is false)', action='store_true',default=False)
parser.add_argument('-v',"--version", help="Print program version", action="version",version='Version '+__version__)
parser.add_argument('-doc',"--documentation", help="Print short program description", action="version",version=__what__)

args = parser.parse_args()

parfile, filename, polelistname, ext_prefix = args.parfile[0], args.infile[0], args.polelist[0], args.ext[0]

#If asked for by user, print sample parameter file
if args.print_parf:
  mgc3_lib.print_sample_parfile()
  sys.exit('Exiting.')

#Parse parameter file
print('Reading Parameter file %s ...' % (parfile))
survey_pars=mgc3_lib.parse_pars(parfile)

#print 'Reading input file %s ...' % (filename)
#obsdata = scipy.genfromtxt(filename,comments='#')
#print 'Input file shape (rows,cols): ', obsdata.shape

print('Reading input file %s ...' % (filename))
obsdata,ugfilename=mgc3_lib.read_inputcat_for_mgc3(filename,pardic=survey_pars)
print('Input file shape (rows,cols): ', obsdata.shape)


#Save input file header
head = myutils.get_header_line(filename)
head[-1]=head[-1]+' IDpole' 

print('Reading pole list file %s ... ' % (polelistname))
polelist=scipy.genfromtxt(polelistname,comments='#',usecols=(0,1,2))
if np.ndim(polelist)==1: polelist=np.array([polelist,])

#Figure the output filename
eind=polelistname.rfind('.'+args.m.lower()) #Find .mgc3,.ngc3,.gc3 extension (the dot is needed to avoid confusion for gc3)
if eind==-1: 
  sys.exit('WARNING: Method extension %s not found in pls filename.\nCheck the consistency of selected method to extract stars with the one used for the pole list provided.\nExiting...' % (args.m.lower()))
if 'peak.pls' in polelistname: pext='peak.pst' 
else: pext='pst'
outfilename='%s%s.%s.%s' % (polelistname[:eind],ext_prefix,args.m.lower(),pext)
print('Printing output file %s ...' % (outfilename))

#Open output file and print the exact same header the input file had, with an extra column
outfile=open(outfilename,'w')
if args.allow_repeats:
 outfile.write('#Allow a star to be printed out multiple times if assoc to >1 pole set to %s\n' % (args.allow_repeats))
outfile.write('#Stars selected according to %s criteria\n' % (args.m))
scipy.savetxt(outfile,head,fmt='%s')

#-------------Solar position and velocity----------------------------
print('Setting solar constants...')
#mycst = my_constants(rsun=8.5,Ugc_hel=10.3,Vgc_hel=220+12.6,Wgc_hel=5.9) --old defaults -changed 25/07/2019
#rsun=8.34, Vc=240 from Reid 2014, U,V,Wsun from Schoenrich & Binney 2010 - as assumed in GaiaCol, Katz et al 2018
#mycst = mgc3_lib.my_constants(rsun=8.,Ugc_hel=11.1,Vgc_hel=240.+12.24,Wgc_hel=7.25) #preferred by Pau - test E
#LM10
#mycst = mgc3_lib.my_constants(rsun=8.,Ugc_hel=9.,Vgc_hel=220.+12.,Wgc_hel=7.) # testLM - used by L&M10
mycst = mgc3_lib.my_constants(rsun=8.,Ugc_hel=11.1,Vgc_hel=235.+12.24,Wgc_hel=7.25) #used in Belokurov2014 -- test D


#The masks associated to each pole will be combined with OR, this way, each star can only be printed once,
#even if associated to more than one pole (unless -allow_repeats is ON)
indep_mask=np.zeros(obsdata[:,0].size,dtype=bool)
indep_pole_ids=-1*np.ones(obsdata[:,0].size)
prev_poleid=1
#Cicle over pole list, one by one
for id_pole,phi_pole,theta_pole in polelist:
  if id_pole!=prev_poleid and args.allow_repeats:
    #Print stored stars 
    print_data=obsdata[indep_mask,:].T
    scipy.savetxt(outfile,np.vstack([print_data,indep_pole_ids[indep_mask]]).T)
    #Reset mask
    indep_mask=np.zeros(obsdata[:,0].size,dtype=bool)
    indep_pole_ids=-1*np.ones(obsdata[:,0].size)
    #Store new poleid
    prev_poleid=id_pole

  mygrid = mgc3_lib.pole_grid(poles=[phi_pole,theta_pole],cst=mycst)
  cat_mask = mygrid.mgc3_allobs_one_pole(obsdata,pars=survey_pars,return_mask=args.m)
  print('   stars associated to pole %s: %d' % (id_pole,np.sum(1*cat_mask)))
  #combine mask with OR
  indep_mask = indep_mask | cat_mask
  #label stars added to the mask in this step with current poleid
  indep_pole_ids[cat_mask] = id_pole

#Print for last pole in the loop
if args.allow_repeats:
  #Print stored stars 
  print_data=obsdata[indep_mask,:].T
  scipy.savetxt(outfile,np.vstack([print_data,indep_pole_ids[indep_mask]]).T)
#Printing is done after finishing the loop if no repetitions are allowed
else:
 print_data=obsdata[indep_mask,:].T
 indep_pole_ids=indep_pole_ids[indep_mask]
 scipy.savetxt(outfile,np.vstack([print_data,indep_pole_ids]).T)

print('Output file written: %s' % (outfilename))


