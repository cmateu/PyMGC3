#!/usr/bin/env python
import numpy as np
import scipy 
import argparse
import os
import sys
import myutils
import gzip

__version__ = '2.0.1'
__docformat__ = "reredtext en"
__what__= sys.argv[0]+": This program combines a list of pole count maps"
#
parser = argparse.ArgumentParser(description='Add mGC3/nGC3/GC3 pole count maps')
parser.add_argument('infile1',metavar='infile1',help='Input pole count map (*.cts file). Gzip supported',nargs=1,action='store')
parser.add_argument('mode1',help='Select mGC3/nGC3/GC3/mGC3hel/GC3hel method for infile1', action='store',choices=['mGC3','nGC3','GC3','mGC3hel','GC3hel'])
parser.add_argument('infile2',metavar='infile2',help='Input pole count map (*.cts file). Gzip supported',nargs=1,action='store')
parser.add_argument('mode2',help='Select mGC3/nGC3/GC3/mGC3hel/GC3hel method for infile2', action='store',choices=['mGC3','nGC3','GC3','mGC3hel','GC3hel'])
parser.add_argument('ofilen',metavar='outfilename',help='Output pole count map name',nargs=1,action='store')
parser.add_argument('-f','--force',help='Force running with existing files, ignoring missing file warnings', action='store_true',default=False)

#---------Parse----------------------------
args = parser.parse_args()

#Read input files
try:
  print 'Reading Input File1:',args.infile1[0]
  #pcm1=scipy.genfromtxt(args.infile1[0])
  if '.gz' in args.infile1[0]: infile1=gzip.open(args.infile1[0],'r')
  else: infile1=open(args.infile1[0],'r')
  pcm1=scipy.genfromtxt(infile1)
except IOError:
  print 'WARNING - File not found: %s' % (args.infile1[0])

try:
  print 'Reading Input File2:',args.infile2[0]
  #pcm2=scipy.genfromtxt(args.infile2[0])
  if '.gz' in args.infile2[0]: infile2=gzip.open(args.infile2[0],'r')
  else: infile2=open(args.infile2[0],'r')
  pcm2=scipy.genfromtxt(infile2)
except IOError:
  print 'WARNING - File not found: %s' % (args.infile1[0])

#Check matching file shapes, exit if mismatch
if np.shape(pcm1)!=np.shape(pcm2):
   print 'Shape1=',np.shape(pcm1)
   print 'Shape2=',np.shape(pcm2)
   sys.exit('WARNING: Input file shapes are inconsistent. Exiting...')

#Initialize all counts as zero
pcm_sum=pcm1.copy()
pcm_sum[:,2:-1]=0

mode1,mode2=args.mode1.lower(),args.mode2.lower()
if 'hel' in mode1 and 'mgc3' in mode1:  counts_col1=4-1
elif 'gc3hel' in mode1: counts_col1=7-1
elif 'mgc3' in mode1: counts_col1=3-1
elif 'ngc3' in mode1: counts_col1=6-1
elif 'gc3'  in mode1: counts_col1=5-1

if 'hel' in mode2 and 'mgc3' in mode2:  counts_col2=4-1
elif 'gc3hel' in mode2: counts_col2=7-1
elif 'mgc3' in mode2: counts_col2=3-1
elif 'ngc3' in mode2: counts_col2=6-1
elif 'gc3'  in mode2: counts_col2=5-1

print 'Print File1 method %s ($%d)' % (mode1,counts_col1+1)
print 'Print File2 method %s ($%d)' % (mode2,counts_col2+1)

print 'Output stored in column $%d (file1 method column)' % (counts_col1+1)

#Combine counts and store them in the column corresp. to the mode selected for input file 1
pcm_sum[:,counts_col1]=pcm1[:,counts_col1]+pcm2[:,counts_col2]
ncountcols=pcm_sum[0,2:-1].size


#Initialize output file
ofilen=args.ofilen[0]
ofile=open(ofilen,'w')
heads=myutils.get_header_line(args.infile1[0])
extra_heads='#%s col = %s_1 + %s_2\n#File1=%s\n#File2=%s\n' % (mode1,mode1,mode2,args.infile1[0],args.infile2[0])
extra_heads=extra_heads+('#'+86*'-')
heads.insert(-1,extra_heads)
for head_line in heads:
  ofile.write(head_line+'\n')

#Printing output file 
countsfmt=ncountcols*'%10d '
print 'Printing output file',ofilen
fmt='%10.3f %10.3f '+countsfmt+'%10.4f' 
scipy.savetxt(ofile,pcm_sum,fmt=fmt)

print 'Done'

