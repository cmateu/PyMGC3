#!/usr/bin/env python
import numpy as np
import scipy 
import argparse
import os
import sys
import myutils

__version__ = '2.0.1'
__docformat__ = "reredtext en"
__what__= sys.argv[0]+": This program combines a list of pole count maps"
#
parser = argparse.ArgumentParser(description='Add mGC3/nGC3/GC3 pole count maps')
parser.add_argument('infilel',metavar='infile_list',help='Input pole count map list (*.cts files)',nargs=1,action='store')
parser.add_argument('ofilen',metavar='outfilename',help='Output pole count map name',nargs=1,action='store')
parser.add_argument('-n','--norm',help='Normalize each PCM *before* adding them together', action='store_true',default=False)
parser.add_argument('-f','--force',help='Force running with existing files, ignoring missing file warnings', action='store_true',default=False)

#---------Parse----------------------------
args = parser.parse_args()

print 'Reading file list',args.infilel
infilelist=scipy.genfromtxt(args.infilel[0],dtype='S',usecols=(0,),unpack=True)
if np.ndim(infilelist)==0: infilelist=np.array([infilelist,])

#Initialize output file
ofilen=args.ofilen[0]
ofile=open(ofilen,'w')
heads=myutils.get_header_line(infilelist[0])
for head_line in heads: 
  ofile.write(head_line+'\n')
#os.system("awk '$0~/#/ {print $0}' %s > %s " % (infilelist[0],ofilen))

#If --norm, print coords for each files's max counts in pls output file
if args.norm:
 omaxfilel=[]
 for ll in ['mgc3','mgc3hel','gc3','ngc3']:
  omaxfilen=ofilen.replace('cts','%s.pls' % (ll))
  omaxfile=open(omaxfilen,'w')
  omaxfile.write('#%5s %10s %10s %10s\n' % ('IDst','phi_max','theta_max','Cnorm'))
  omaxfilel.append(omaxfile)

for n in range(len(infilelist)):
   
   infile=infilelist[n]

   print 'Reading file %d of %d (%s)' % (n+1,len(infilelist),infile)
   try:
    pcm=scipy.genfromtxt(infile)
   except IOError:
    print 'WARNING - File not found: %s' % (infile)
    if args.force: 
      print 'Continue...'
      continue
    else: sys.exit('Exiting (to force skipping missing files, use -f option)')

   ncountcols=pcm[0,2:-1].size
   nfinalccol=2+ncountcols

   #Initialize normalization factor (off by default)
   fnorm_vec=1.
   if args.norm: 
     fnorm_vec=1./np.max(pcm[:,2:nfinalccol],axis=0).astype(float)
     fmax_ind=np.argmax(pcm[:,2:nfinalccol],axis=0)
     fnorm_vec[fnorm_vec==0.]=1.  #If the max of any column is zero, change multiplication factor to 1.
     for kk in range(fmax_ind.size):
      try: iidst=infile[infile.find('id')+2:infile.find('id')+2+3]
      except: iidst=n
      omaxfilel[kk].write('%6s ' % (iidst)) 
      omaxfilel[kk].write('%10.3f %10.3f %10d\n' % tuple(pcm[fmax_ind[kk],[0,1,kk+2]])) 
       
   if n==0:
     pcm_sum=pcm  #Initialize matrix with data for the first file
     pcm_sum[:,2:nfinalccol]=pcm_sum[:,2:nfinalccol]*fnorm_vec   #Normalize each of the counts columns
     pcm_shape=np.shape(pcm)   #Save first file's shape to check consistency with the rest of the files
   else:
     #Check shape
     if pcm_shape!=np.shape(pcm):
        print 'Shape nf= 1',pcm_shape
        print 'Shape nf=',n+1,np.shape(pcm)
        sys.exit('WARNING: Input file shapes are inconsistent. Exiting...')
     #If not first file, add the columns corresponding to pole counts. Leave the rest as in the first file
     pcm_sum[:,2:nfinalccol]=pcm_sum[:,2:nfinalccol]+ (pcm[:,2:nfinalccol]*fnorm_vec) 
  

#Printing output file 
countsffmt=ncountcols*'%10.4g '
countsfmt=ncountcols*'%10d '
print 'Printing output file',ofilen
if args.norm:
  fmt='%10.3f %10.3f '+countsffmt+'%10.4f' 
  scipy.savetxt(ofile,pcm_sum,fmt=fmt)
else:
  fmt='%10.3f %10.3f '+countsfmt+'%10.4f' 
  scipy.savetxt(ofile,pcm_sum,fmt=fmt)

print 'Done'

