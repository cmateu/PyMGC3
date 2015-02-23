#!/usr/bin/env python
import numpy as np
import scipy 
import argparse
import os
import sys

__version__ = '2.0.1'
__docformat__ = "reredtext en"
__what__= sys.argv[0]+": This program combines a list of pole count maps"
#
parser = argparse.ArgumentParser(description='Add mGC3/nGC3/GC3 pole count maps')
parser.add_argument('infilel',metavar='infile_list',help='Input pole count map list (*.cts files)',nargs=1,action='store')
parser.add_argument('ofilen',metavar='outfilename',help='Output pole count map name',nargs=1,action='store')

#---------Parse----------------------------
args = parser.parse_args()

print 'Reading file list',args.infilel
infilelist=scipy.genfromtxt(args.infilel[0],dtype='S')

#Initialize output file
ofilen=args.ofilen[0]
os.system("awk '$0~/#/ {print $0}' %s > %s " % (infilelist[0],ofilen))
ofile=open(ofilen,'a')

for n in range(len(infilelist)):
   
   infile=infilelist[n]

   print 'Reading file %d of %d (%s)' % (n+1,len(infilelist),infile)
   pcm=scipy.genfromtxt(infile)

   if n==0:
     pcm_sum=pcm  #Initialize matrix with data for the first file
     pcm_shape=np.shape(pcm)   #Save first file's shape to check consistency with the rest of the files
   else:
     #Check shape
     if pcm_shape!=np.shape(pcm):
        print 'Shape nf= 1',pcm_shape
        print 'Shape nf=',n+1,np.shape(pcm)
        sys.exit('WARNING: Input file shapes are inconsistent. Exiting...')
     #If not first file, sum the columns corresponding to pole counts. Leave the rest as in the first file
     pcm_sum[:,2:6]=pcm_sum[:,2:6]+pcm[:,2:6]   

#Printing output file 
print 'Printing output file',ofilen
scipy.savetxt(ofile,pcm_sum,fmt='%10.3f %10.3f %10d %10d %10d %10d %10.4f')
print 'Done'

