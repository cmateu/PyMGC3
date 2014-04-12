#!/usr/bin/env python

import Image
import sys
from matplotlib.mlab import griddata
import pylab as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.image import NonUniformImage
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage import maximum_position
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
from scipy.ndimage import median_filter
from scipy.ndimage import label
from scipy.ndimage import center_of_mass
import scipy
import scipy.interpolate
import argparse

def get_map4plot(lpole,bpole,npole):
 lp,bp=np.radians(lpole),np.radians(bpole)
 mask=lp>np.pi
 lp[mask]=lp[mask]-2*np.pi
 #duplicate to add southern hemisphere
 lp1=lp+np.pi
 mask=lp1>np.pi
 lp1[mask]=lp1[mask]-2*np.pi
 lp2=np.append(lp,lp1[bp>0.])
 bp2=np.append(bp,-bp[bp>0.])
 np2=np.append(npole,npole[bp>0.])

 return (lp2,bp2,np2)

def detect_peaks(image,nsigma=3.,size=10,nsigclip=3.):
    """
    Takes an image and detect the peaks using the local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """

    # define an 8-connected neighborhood
    neighborhood = generate_binary_structure(2,2)

    #apply the local maximum filter; all pixel of maximal value in their neighborhood are set to 1
    #local_max = maximum_filter(image, footprint=neighborhood)==image
    local_max = maximum_filter(image, size=10, mode='wrap')==image

    mean=np.mean(image)
    std=np.std(image)
    cmean=np.mean(image[np.abs(image-mean)<=nsigclip*std])
    sigma=np.std(image[np.abs(image-mean)<=nsigclip*std])
    #print 'mean',cmean,'sig',sigma
    background = (image<(cmean+nsigma*sigma))  #image==0

    #a little technicality: we must erode the background in order to 
    #successfully subtract it form local_max, otherwise a line will 
    #appear along the background border (artifact of the local maximum filter)
    eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)

    #we obtain the final mask, containing only peaks, 
    #by removing the background from the local_max mask
    detected_peaks = local_max - eroded_background

    #Extract positions for identified maxima
    im_labels,im_num = label(detected_peaks)
    print 'Detections:', im_num
    if im_num==0:
      return []
    pos = center_of_mass(detected_peaks,im_labels,range(im_num))


    return (detected_peaks,np.array(pos))

def interpolate_peak_location(peaks,xgrid,ygrid):

 #Parse peaks matrix
 xp_ind=peaks[1:,1]
 yp_ind=peaks[1:,0]

 #Create interpolators
 xinds=range(lgrid.size)
 yinds=range(bgrid.size)
 getx_from_ind = scipy.interpolate.interp1d(xinds,xgrid,kind='linear')
 gety_from_ind = scipy.interpolate.interp1d(yinds,ygrid,kind='linear')

 #Evaluate
 xp=getx_from_ind(xp_ind)
 yp=gety_from_ind(yp_ind)

 return (xp,yp)

#---------------------------------------------------------

#Parse command line arguments

__version__ = '1.0'
__docformat__ = "reredtext en"
__what__= sys.argv[0]+": This program detects peaks in an input mgc3 pole count map, prints them out and outputs a pole count map with the detections indicated"
#
parser = argparse.ArgumentParser()
#
parser.add_argument('infile',metavar='infile.mgc3.cts',help='Input file containing mgc3-counts map',nargs=1,action='store')
parser.add_argument('-ngrid',help='Number of pixels (in 1-D) for interpolated image. Default=300', nargs='?',default=300, type=int)
parser.add_argument('-msize',help='Neighbourhood size (in pixels) for median filter. Default=60', nargs='?',default=60,type=int)
parser.add_argument('-nsigma',help='Nsigma threshold for peak detection. Default=3', nargs='?', default=3,type=float)
parser.add_argument('-psize',help='Rough estimate of peak size or extent (in pixels). Default=10', nargs='?',default=10,type=int)
parser.add_argument('-gc3',help='Use gc3 instead of mgc3 counts', action='store_true')
parser.add_argument('-hel',help='Use heliocentric instead of galactocentric counts', action='store_true')
parser.add_argument('-farea',help='For non-all sky surveys, correct counts by area fraction', action='store_true')
parser.add_argument('-s', help='Do not display output pole count map',action='store_false')
parser.add_argument('-sa','--show_all',help='Display all maps used (i.e. smoothed and subtracted)',action='store_true')
parser.add_argument('-v',"--version", help="Print program version", action="version",version='Version '+__version__)
parser.add_argument('-doc',"--documentation", help="Print short program description", action="version",version=__what__)

args = parser.parse_args()

#-------------------------------------------------------------------------------------------------------

#Output filename
infile=args.infile[0]
outfile=infile.replace('.mgc3.cts','')+'.mgc3.pls'

#Choose the appropriate column depending on the input counts chosen by user
if not args.gc3 and not args.hel: nc=2
if not args.gc3 and args.hel: nc=3
if args.gc3 and not args.hel: nc=4
if args.gc3 and args.hel: sys.exit('GC3-hel counts not available')

lpole,bpole,npole,farea=scipy.genfromtxt(infile,unpack=True,usecols=(0,1,nc,5))

if args.farea: 
  print 'Correcting for survey area coverage...'
  npole=npole/farea #Default is false

d2r=np.pi/180.
lp,bp,npole=get_map4plot(lpole,bpole,npole)

#Make image
lgrid=np.linspace(-np.pi,+np.pi,args.ngrid)
bgrid=np.linspace(-np.pi/2.,+np.pi/2.,args.ngrid)
im=griddata(lp,bp,npole,lgrid,bgrid,interp='linear')

#Image props
cmap=cm.jet
cmap=cm.spectral
extent=[-np.pi,+np.pi,-np.pi/2.,+np.pi/2.]

#Smoothed image
smooth=median_filter(im,size=args.msize)
#Unsharp-masked image
subim=im-smooth

fig1=plt.figure(1,figsize=(13,8))
ax1=fig1.add_subplot(221) #,projection='aitoff')
sc=ax1.imshow(im,interpolation='bilinear',aspect='auto',origin='lower',cmap=cmap,extent=extent)
plt.colorbar(sc,pad=0)

#Detect peaks
detected_peaks_im, peaks = detect_peaks(im,nsigma=args.nsigma,nsigclip=3.,size=args.psize)
lpeaks,bpeaks = interpolate_peak_location(peaks,lgrid,bgrid)
lpeaks, bpeaks = lpeaks[bpeaks>=0.], bpeaks[bpeaks>=0.]

#Smoothed image
ax2 = fig1.add_subplot(222)
sc2=ax2.imshow(smooth,cmap=cmap,origin='lower',aspect='auto',extent=extent)
plt.colorbar(sc2,pad=0)

#Unsharp masked image
ax3 = fig1.add_subplot(223)
sc3=ax3.imshow(subim,cmap=cmap,origin='lower',aspect='auto',extent=extent)
plt.colorbar(sc3,pad=0)

ax4 = fig1.add_subplot(224)
sc4=ax4.imshow(subim,cmap=cmap,origin='lower',aspect='auto',extent=extent)
ax4.plot(lpeaks,bpeaks,'w^',ms=8,mfc='None',mec='white')
ax4.set_xlim(extent[0],extent[1])
ax4.set_ylim(extent[2],extent[3])

#Show auxiliary plot
if args.show_all: plt.show()
plt.close()

#Pole count map -----------------------------------------------
fig2=plt.figure(2,figsize=(12,6))
ax = fig2.add_subplot(111,projection='aitoff')
sc=ax.scatter(lp,bp,c=npole,edgecolors='none',s=(40+15./np.cos(bp)),cmap=cm.spectral)
ax.plot(lpeaks,bpeaks,'w^',ms=7) #,mfc='None',mec='white',lw=3)
plt.savefig(outfile.replace('pls','png'))
if args.s: plt.show()

#Print pole coordinates in output file (only for northern hemisphere)
lpeaks[lpeaks<0.]=lpeaks[lpeaks<0.]+2*np.pi
ids=range(lpeaks.size)
outf=open(outfile,'w')
scipy.savetxt(outf,np.array(['#ID phi_pole theta_pole',]),fmt='%s')
scipy.savetxt(outf,np.array([ids,lpeaks/d2r,bpeaks/d2r]).T,fmt='%3d %8.4f %8.4f')
