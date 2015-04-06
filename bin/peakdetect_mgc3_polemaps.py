#!/usr/bin/env python
import pylab as plt
import scipy
import numpy as np
from mpl_toolkits.basemap import Basemap
import sys
import argparse
import scipy.ndimage
import pyfits
import os

class xypix_converter:
  
  def __init__(self,m,npix=256,rangex=(0.,1.),rangey=(0.,1.)):

    self.xo,self.xf=rangex
    self.yo,self.yf=rangey
    self.npix=npix

    self.dx=(self.xf-self.xo)/np.float(self.npix)
    self.dy=(self.yf-self.yo)/np.float(self.npix)

    self.basemap_trans=m

  def get_phys_from_pix(self,xpix,ypix):

    xc=(xpix)*self.dx + self.xo
    yc=(ypix)*self.dy + self.yo

    #Convert cartesian to spherical coords  phic,thetac=m(xc,yc,inverse=True)
    phic,thetac=self.basemap_trans(xc,yc,inverse=True)
    mask=phic<0.
    phic[mask]=phic[mask]+360.

    return (xc,yc,phic,thetac)

#----------------------------------------------

__version__ = '2.0.1'
__docformat__ = "reredtext en"
__what__= sys.argv[0]+": This program detects peaks in pole count maps using the Fellwalker algorithm (starlink implementation)"
#
parser = argparse.ArgumentParser(description='Detect peaks in mGC3/nGC3/GC3 pole count maps')
parser.add_argument('infile',metavar='infile',help='Input file containing pole count maps (*.cts file)',nargs=1,action='store')
parser.add_argument("-l", "--llist", action="store_true",help='Take infile as list of mgc3.cts files')
parser.add_argument('-m',help='Plot mGC3/nGC3/GC3 pole count map. Default is mGC3', action='store',default='mGC3',choices=['mGC3','nGC3','GC3'])
parser.add_argument('-f','--fig',help='Output plot type png/eps. Default is png', action='store',default='png',choices=['png','eps','pdf'])
parser.add_argument('-ext',metavar='outfile_ext',help='Output suffix [optional]. If given output will be infile.outfile_ext.mgc3.pst',action='store',default=['',],nargs=1)
parser.add_argument('-log',help='Plot detected peaks in log-count map', action='store_true',default=False)
parser.add_argument('-labels',help='Plot peak ID labels', action='store_true',default=False)
parser.add_argument('-title',help='Plot title', action='store',default=None)
parser.add_argument('-lon0',help='Longitude for Y-axis. Default is 0.', action='store',default=0.,type=np.float)
parser.add_argument('-lat0',help='Bounding latitude for plot. Default is 90.', action='store',default=0.,type=np.float)
parser.add_argument('-vmin',help='Min counts for color-scale. Default is min(cts)', action='store',default=None,type=np.float)
parser.add_argument('-vmax',help='Max counts for color-scale. Default is max(cts)', action='store',default=None,type=np.float)
parser.add_argument('-dlat',help='Spacing between parallels. Default is 30.', action='store',default=20.,type=np.float)
parser.add_argument('-dlon',help='Spacing between meridians. Default is 30.', action='store',default=30.,type=np.float)
parser.add_argument('-ms',help='Marker size. Default: 50 for npaeqd.', action='store',default=50,type=np.float)
parser.add_argument('-r','--raw',help='Plot raw grid pole-count map instead of contour map.', action='store_true',default=False)
parser.add_argument('-t','--twohemispheres',help='Plot both hemispheres in pole-count map.', action='store_true',default=False)
parser.add_argument('-s','--show',help='Show plot in window. Default is False', action='store_true',default=False)
parser.add_argument('-sc','--saveclumps',help='Plot and save poles associated to each peak.', action='store_true',default=False)
parser.add_argument('-mj','--maxjump',help='Fellwalker MaxJump param, neighbourhood radius to search for +gradient. Default 20.', action='store',default=20,type=np.float)
parser.add_argument('-al','--alpha',help='Clump transparency. Default 0.3', action='store',default=0.3,type=np.float)
peakargs = parser.add_mutually_exclusive_group()
peakargs.add_argument('-fr','--frms',help='If set, min peak height is frms*RMS', action='store',type=np.float)
peakargs.add_argument('-ff','--ffrac',help='Default option. Min peak height is fmax*max_pole_counts. Default fmax=0.6', action='store',default=0.6,type=np.float)
parser.add_argument('-fx','--fwxm',help='Store pixels with cts>fwxm*maxpeak. Default is 0.5 (=FWHM)', action='store',default=0.5,type=np.float)

parser.add_argument('-U','--unsharp',help='Detect peaks in unsharp-masked image. Default is False', action='store_true',default=False)
parser.add_argument('-ns','--nsigma',help='If -U is set, set N-sigma threshold for detections. Default is 3', action='store',type=np.float,default=3.)


#---------Parse----------------------------
args = parser.parse_args()

#Parse inputs

if not args.llist:
 print 'Reading file: %s' % (args.infile)
 file_list=[args.infile[0],]
else:
 print 'Reading input files from list file: ', args.infile[0]
 file_list=scipy.genfromtxt(args.infile[0],dtype='S')
 if np.ndim(file_list)==0: file_list=array([file_list,])

#Mode-------------------------------------------
mode=args.m.lower()
mode_ori=args.m
print 'Pole counts plotted: ', mode_ori
if 'mgc3' in mode:   counts_col=3-1
elif 'ngc3' in mode: counts_col=6-1
elif 'gc3'  in mode: counts_col=5-1

#Parse raw/contour mode-------------------------
if args.raw: 
  pmode='r'
  print 'Plotting raw pole-count map'
else: 
  pmode='c'
  print 'Plotting contour pole-count map'

if args.twohemispheres:
  print 'Plotting both hemispheres in pole-count map'
else: print 'Plotting one hemisphere in pole-count map.'

proj='npaeqd'
print 'Plotting using projection:', proj

ori='vertical'
ori='horizontal'
ni=0

colormap=plt.cm.jet
colormap=plt.cm.gray
for infilen in file_list:

  phio,thetao,pole_ctso=pdat=scipy.genfromtxt(infilen,comments='#',usecols=(0,1,counts_col),unpack=True)
  #If log-flag is set, do everything with log(counts)
  if args.log: 
     print 'Do peak detection on linear scale but show plot on log scale'
     pmode=pmode+'l'

  #Default title----------------------------------
  args.title=infilen

  #Output figure and file names
  figname_root=infilen.replace('.mgc3.cts',args.ext[0])  #works well if args.ext is empty
  figname='%s.%s.%s.%s.%s' % (figname_root,mode,proj[:3],pmode,args.fig)
  clumpfname='%s.%s.peak.pls' % (figname_root,mode)
  clumppixfname='%s.%s.pls' % (figname_root,mode)
  print 'Output filename:', figname

  phi2=np.append(phio,(phio+180.) % 360.)
  theta2=np.append(thetao,-thetao)
  pole_cts2=np.append(pole_ctso,pole_ctso) 

  if args.twohemispheres or args.unsharp:
   phi,theta,pole_cts=phi2,theta2,pole_cts2
  else: 
   phi,theta,pole_cts=phio,thetao,pole_ctso 
  #This is done this way so that for contour plotting, the twohemispheres are used for the interpolation
  #into a regular grid, but only one is shown 

  #----------------Pole count map------------------------------
  mer_grid=[0.,360.,args.dlon]
  par_grid=[-args.dlat,+90.,args.dlat]

  #PCM map projection basics
  fig=plt.figure(1,figsize=(8,8))
  dw=0.8
  wo=(1.-dw)/2.
  wyo=0.75*wo 
  wyo=0.05*wo 
  fig.subplots_adjust(left=wo,right=dw+wo,top=dw+wo,bottom=wyo)
  nrow,ncol,nplot=1,1,1
  l0=args.lon0
  proj_dict={'boundinglat':args.lat0,'resolution':'l','lon_0':l0}
  ms=args.ms

  ax=fig.add_subplot(nrow,ncol,nplot)
  m = Basemap(projection=proj,ax=ax,**proj_dict)
  m.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]),color='lightgrey',lw=2.)
  m.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]),color='lightgrey',lw=2.)
  m.drawmapboundary()
  x,y=m(phi,theta)

  #------------Grid-data for contour plotting---------------
  npix=400
  clevels=30
  xo,xf=np.min(x),np.max(x)
  yo,yf=np.min(x),np.max(x)
  xi = np.linspace(xo,xf,npix)
  yi = np.linspace(yo,yf,npix)
  zi = plt.griddata(x,y,pole_cts,xi,yi) #,'nn')

  lmax=np.floor(np.log10(np.max(pole_cts)))
  if 'r' in pmode: 
     if args.log: c=m.scatter(x,y,c=np.log10(pole_cts),edgecolor='none',s=ms,cmap=colormap,vmin=args.vmin,vmax=args.vmax)
     else:        
        if args.vmin is not None: vmin=args.vmin/10**lmax
        else: vmin=args.vmin
        if args.vmax is not None: vmax=args.vmax/10**lmax
        else: vmax=args.vmax
        c=m.scatter(x,y,c=pole_cts/10**lmax, edgecolor='none',s=ms,cmap=colormap,vmin=vmin,vmax=vmax)
  else: 
   if args.log: c=m.contourf(xi,yi,np.log10(zi),clevels,cmap=colormap,vmin=args.vmin,vmax=args.vmax)
   else:
     #Use vmin and vmax for plotting only
     if args.vmin is not None: vmin=args.vmin
     else: vmin=np.min(zi)
     if args.vmax is not None: vmax=args.vmax
     else: vmax=np.max(zi)
     zii=zi
     zii[(zii<vmin)]=vmin
     zii[(zii>vmax)]=vmax
     lmax=np.floor(np.log10(vmax))
     c=m.contourf(xi,yi,zii/10**lmax, clevels,cmap=colormap)

  #Plot colorbar
  cax0=plt.gca().get_position()
  cax=plt.axes([cax0.x0,cax0.y0+dw+0.05,dw,0.02])
  cb=plt.colorbar(c,cax=cax,orientation='horizontal',format='%4.1f')
  cax.xaxis.set_ticks_position('top')
  
  #Labels and such
  if lmax>0: factorl='$\\times 10^{%d}$ ' % (lmax)
  else: factorl=''
  if args.log:
     cax.set_xlabel('%s log-pole-counts ($\log_{10}$-stars/pole)' % (mode_ori))
  else:
     cax.set_xlabel('%s pole-counts (%sstars/pole)' % (mode_ori,factorl))
  cax.xaxis.set_label_position('top') 

  if args.title:
    ax.text(0.5,1.14,args.title,transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontsize=16)

  #-------------------------Unsharp masking-----------------------------------------------------
  if args.unsharp:
    #Compute median filtered image
    nsm=60 #neighbourhood for median computation
    zi_smooth=scipy.ndimage.median_filter(zi,size=(nsm,nsm),mode='wrap')
    zi_sharp=zi-zi_smooth
    zi_sigma=np.sqrt(zi_smooth)  #assuming counts in the PCM follow a Poisson distribution
    zi_sharp_Nsig=zi_sharp/zi_sigma   #express subtracted image in nsigma-units
    #----plot unsharp masked and smoothed images
    fig2=plt.figure(2,figsize=(16,6))
    #----------------smoothed image---------------------------------
    axs=fig2.add_subplot(1,3,1)
    ms = Basemap(projection=proj,ax=axs,**proj_dict)
    ms.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]),color='lightgrey',lw=2.)
    ms.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]),color='lightgrey',lw=2.)
    ms.drawmapboundary()
    c1=ms.contourf(xi,yi,zi_smooth, clevels,cmap=colormap)
    cb=plt.colorbar(c1,ax=axs,orientation='horizontal',format='%d',pad=0,aspect=30)
    #-------------subtracted image-----------------------------------------
    axu=fig2.add_subplot(1,3,2)
    mu = Basemap(projection=proj,ax=axu,**proj_dict)
    mu.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]),color='lightgrey',lw=2.)
    mu.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]),color='lightgrey',lw=2.)
    mu.drawmapboundary()
    c3=mu.contourf(xi,yi,np.log10(zi_sharp),clevels,cmap=colormap,vmin=0.)
    cb=plt.colorbar(c3,ax=axu,orientation='horizontal',format='%4.1f',pad=0,aspect=30)
    #----------------subtracted image in Nsgima units---------------------------------
    axn=fig2.add_subplot(1,3,3)
    mn = Basemap(projection=proj,ax=axn,**proj_dict)
    mn.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]),color='lightgrey',lw=2.)
    mn.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]),color='lightgrey',lw=2.)
    mn.drawmapboundary()
    colormap_nsig=plt.cm.spectral
    #colormap_nsig.set_over('firebrick')
    zi_sharp_cut=zi_sharp_Nsig.copy()
    sigmax=12.  #Maximum Nsigma for contour and color display
    zi_sharp_cut[zi_sharp_cut>=sigmax]=sigmax
    c2=mn.contourf(xi,yi,(zi_sharp_cut), np.arange(0.,sigmax+1.,1.),cmap=colormap_nsig)
    cb=plt.colorbar(c2,ax=axn,orientation='horizontal',format='%4.1f',pad=0,aspect=30,extend='max')
    #---------
    fig2.tight_layout()
    print 'here----'
    usharp_figname='%s.%s.%s.%s.usharp.%s' % (figname_root,mode,proj[:3],pmode,args.fig)
    fig2.savefig(usharp_figname)

  #--------------------Manage formats to run the Fellwaker peak detection algorithm---------------------------------------------
  if args.unsharp: zi=zi_sharp_Nsig
 
  hdu = pyfits.PrimaryHDU()
  zi.data[zi is np.nan]=-1
  hdu.data=zi.data
  hdu.writeto('_zp.fits', clobber=True)

  #convert fits to ndf
  starlink_path='/star-2014A/bin/'
  print 'Converting fits to NDF'
  os.system('%s/convert/fits2ndf _zp.fits _zpndf' % (starlink_path)) 
  print '#--------stats-------------------------'
  print 'Mean,Median:',np.mean(zi), np.median(zi) 
  print 'Sqrt(mean), sqrt(median):',np.sqrt(np.mean(zi)), np.sqrt(np.median(zi))
  print '25,50,75 pcntiles:',np.percentile(zi,25), np.percentile(zi,50),np.percentile(zi,75)
  print '#----------------------------------------'
  if args.unsharp:
     rms=1.  #If -U, zi is by always in N-sigma units by default            
     minheight=args.nsigma
  else:
   rms=np.std(zi)
   if args.frms is not None:
     minheight=args.frms*rms
   else: 
     minheight=args.ffrac*np.max(zi)
  print 'Finding clumps with Fellwalker'
  os.system('%s/cupid/findclumps in=_zpndf.sdf out=_zp_cmask method=fellwalker outcat=_zp_clumps rms=%f config="fellwalker.maxjump=%.0f" ' % (starlink_path,rms,args.maxjump)) 

  #------------------------------Deal with identified peaks (peak centroids, etc)---------------------------
  #Read-in output table with identified clumps
  clumpdat=pyfits.open('_zp_clumps.FIT')[1].data                    
  xpix,ypix,cheight=clumpdat.field('Peak1'),clumpdat.field('Peak2'),clumpdat.field('Peak')
  xpix_c,ypix_c,dx_pix,dy_pix=clumpdat.field('Cen1'),clumpdat.field('Cen2'),clumpdat.field('Size1'),clumpdat.field('Size2')
  pid=np.arange(xpix.size) + 1  #Peak IDs, must start from 1

  #Keep only peaks with counts > minheight and with sizes>=1 (clumps have to be at least 1 pixel in size)
  mask=(cheight>=minheight) & (dx_pix>=1.) & (dy_pix>=1.) #& (thetapeak>=0.)
  if not mask.any(): 
    print 'No peaks above minheight or with size>=1pix found'
    continue
  xpix,ypix,cheight=xpix[mask],ypix[mask],cheight[mask]
  xpix_c,ypix_c,dx_pix,dy_pix,pid=xpix_c[mask],ypix_c[mask],dx_pix[mask],dy_pix[mask],pid[mask]
  newid=np.arange(xpix.size) + 1  #Rename so IDs will be consecutive numbers starting from 1

  #Convert pix coords to physical units 
  pix_convert=xypix_converter(m,npix=npix,rangex=(xo,xf),rangey=(yo,yf))
  #Convert peak coords
  xpeak,ypeak,phipeak,thetapeak=pix_convert.get_phys_from_pix(xpix,ypix)

  #---------------Get All Pixels associated to each clump-------------------------------------
  #Convert clump pixel mask to fits
  print 'Converting pixel map from NDF to fits'
  os.system('%s/convert/ndf2fits _zp_cmask _zp_cmask.fits' % (starlink_path))

  #Read in clump-mask fits image
  cmask_dat=pyfits.getdata('_zp_cmask.fits',0)
  inds=np.arange(npix)+1
  xpix_2d,ypix_2d=np.meshgrid(inds,inds)
  xinds,yinds=xpix_2d.flatten(),ypix_2d.flatten()
  pcts_1d=zi.flatten()
  cmask_1d=cmask_dat.flatten() #cmask goes with pid, i.e. pid=cmask_1d.unique()
  xcmask,ycmask,phicmask,thetacmask = pix_convert.get_phys_from_pix(xinds,yinds)

  #----------Deal with clumps that go from one hemisphere to the other-------------------------------
  #Loop over clump IDs, if centroid is in the N, save all its pixels
  #Peak data
  u_phipeak,u_thetapeak,u_xpix,u_ypix=np.array([]),np.array([]),np.array([]),np.array([])
  u_dx_pix,u_dy_pix,u_pid=np.array([]),np.array([]),np.array([])
  #Pixel data
  u_xcmask,u_ycmask=np.array([]),np.array([])
  u_phicmask,u_thetacmask=np.array([]),np.array([])
  u_cmask_1d,u_pcts_1d=np.array([]),np.array([])
  for ii in range(thetapeak.size):
   if thetapeak[ii]>=0: #North peaks
     #Save peak data
     u_phipeak,u_thetapeak=np.append(u_phipeak,phipeak[ii]),np.append(u_thetapeak,thetapeak[ii])
     u_xpix,u_ypix=np.append(u_xpix,xpix[ii]),np.append(u_ypix,ypix[ii])   
     u_dx_pix,u_dy_pix=np.append(u_dx_pix,dx_pix[ii]),np.append(u_dy_pix,dy_pix[ii])   
     u_pid=np.append(u_pid,pid[ii])
     #Save pixel-peak data
     peakmask=cmask_1d==pid[ii]    
     u_xcmask,u_ycmask=np.append(u_xcmask,xcmask[peakmask]),np.append(u_ycmask,ycmask[peakmask])
     u_phicmask,u_thetacmask=np.append(u_phicmask,phicmask[peakmask]),np.append(u_thetacmask,thetacmask[peakmask])
     u_cmask_1d,u_pcts_1d=np.append(u_cmask_1d,cmask_1d[peakmask]),np.append(u_pcts_1d,pcts_1d[peakmask])
   
  from astropy.coordinates import SkyCoord
  from astropy import units as aunits

  #Loop over clump with centroids in the S. If all pixels are in the S, dont save, if some in the N, flip and match
  for ii in range(thetapeak.size):
   if thetapeak[ii]<0: #South peaks
     #Check whether all pixels for a given peak are in the south
     peakmask=(cmask_1d==pid[ii]) 
     smask=thetacmask[peakmask]>=0
     if smask.any():
        print 'Clump ID=%d has theta>0 pixels' % (pid[ii])
        #Flip coords 
        phic,thetac= (phicmask[peakmask]+180.) % 360.,-thetacmask[peakmask]
        cmaskc=cmask_1d[peakmask]
        catN = SkyCoord(u_phicmask*aunits.degree,u_thetacmask*aunits.degree,frame='icrs') 
        catS = SkyCoord(phic*aunits.degree,thetac*aunits.degree,frame='icrs') 
        Nindex,dist2d,dist3d = catS.match_to_catalog_sky(catN) #Returns catN indices for objs matching catS
        mask_tol=dist2d.deg<1. #Keep only matches within less than 1deg
        cmaskc[mask_tol]=u_cmask_1d[Nindex[mask_tol]]
        #Put everything back together 
        u_xcmask,u_ycmask=np.append(u_xcmask,xcmask[peakmask][mask_tol]),np.append(u_ycmask,ycmask[peakmask][mask_tol])
        u_phicmask,u_thetacmask=np.append(u_phicmask,phicmask[peakmask][mask_tol]),np.append(u_thetacmask,thetacmask[peakmask][mask_tol])
        print xcmask[peakmask][mask_tol].size, phicmask[peakmask][mask_tol].size
        print cmaskc.size, pcts_1d[peakmask][mask_tol].size
        u_cmask_1d,u_pcts_1d=np.append(u_cmask_1d,cmaskc[mask_tol]),np.append(u_pcts_1d,pcts_1d[peakmask][mask_tol])

#  #Match detections in the S to those in the N (N is the reference one, S is affected by projection for theta<-30deg)
#  catN = SkyCoord(phi_n*aunits.degree,theta_n*aunits.degree,frame='icrs') 
#  catS = SkyCoord(phi_s*aunits.degree,theta_s*aunits.degree,frame='icrs') 
#  Nindex,dist2d,dist3d = catS.match_to_catalog_sky(catN) #Returns catN indices for objs matching catS
 
     else:
       continue

  #Just for tests
  xcmask,ycmask=u_xcmask,u_ycmask
  phicmask,thetacmask=u_phicmask,u_thetacmask
  cmask_1d,pcts_1d=u_cmask_1d,u_pcts_1d 

#  #Separate data into North and South hemispheres
#  mask=(thetacmask>=0.) & (thetacmask<=90.) & (pcts_1d>=minheight)
#  phi_n,theta_n,cmask_n=phicmask[mask],thetacmask[mask],cmask_1d[mask]
#  xcmask_n,ycmask_n,pcts_1d_n=xcmask[mask],ycmask[mask],pcts_1d[mask]
#  mask=(thetacmask<0.) & (thetacmask>=-90) & (pcts_1d>=minheight)
#  phi_s,theta_s,cmask_s=phicmask[mask],thetacmask[mask],cmask_1d[mask]
#  xcmask_s,ycmask_s,pcts_1d_s=xcmask[mask],ycmask[mask],pcts_1d[mask]
#  #Flip southern hemisphere to match with northern hemisphere data
#  phio_s,thetao_s=phi_s,theta_s #keep unflipped coords for later
#  phi_s,theta_s=phi_s + 180. ,-theta_s
#  phi_s[phi_s>=360.]=phi_s[phi_s>=360.]-360.
#  #Match detections in the S to those in the N (N is the reference one, S is affected by projection for theta<-30deg)
#  from astropy.coordinates import SkyCoord
#  from astropy import units as aunits
#  catN = SkyCoord(phi_n*aunits.degree,theta_n*aunits.degree,frame='icrs') 
#  catS = SkyCoord(phi_s*aunits.degree,theta_s*aunits.degree,frame='icrs') 
#  Nindex,dist2d,dist3d = catS.match_to_catalog_sky(catN) #Returns catN indices for objs matching catS
#  print '---------->',aunits.degree
#  print dist2d.deg
#  #Set S-pixel clumpIDs to those of their N-counterparts
#  mask_tol=dist2d.deg<1. #Keep only matches within less than 1deg
#  cmask_s[mask_tol]=cmask_n[Nindex[mask_tol]]
#  #Put everything back together 
#  xcmask,ycmask=np.append(xcmask_n,xcmask_s[mask_tol]),np.append(ycmask_n,ycmask_s[mask_tol])
#  phicmask,thetacmask=np.append(phi_n,phio_s[mask_tol]),np.append(theta_n,thetao_s[mask_tol])
#  cmask_1d,pcts_1d=np.append(cmask_n,cmask_s[mask_tol]),np.append(pcts_1d_n,pcts_1d_s[mask_tol])
  
  #----------------------Printing Ouput peak-file Header and Params on Screen----------------------------------
  #Open output file to store clump coords
  clumpfile=open(clumpfname,'w')
  #Save some param data and print on screen and as file header
  header_info='#----------------------------------------------------------------------\n'
  header_info=header_info+'# Peak-detection algorithm: Starlink Fellwalker (Berry 2014)\n'
  if args.unsharp: header_info=header_info+'# Params: RMS=%.1f (N-sigma units)\n' % (rms)
  else:            header_info=header_info+'# Params: RMS=%.1f (=sqrt(mean(cts))\n' % (rms)
  header_info=header_info+'#         FellWalker.MaxJump=%.0f\n' % (args.maxjump)
  if args.unsharp:
    header_info=header_info+'#         MinHeight=%-4.1f (Nsigma-units)\n' % (minheight)
  else:
   if args.frms is not None: 
     header_info=header_info+'#         MinHeight=%d (=%.1f*RMS)\n' % (minheight,args.frms)
   else: 
    header_info=header_info+'#         MinHeight=%d (=%.2f*max_cts)\n' % (minheight,args.ffrac)
  header_info=header_info+'#         FWXM=%.2f (stored pixels associated to each clump)\n' % (args.fwxm)
  header_info=header_info+'#----------------------------------------------------------------------\n'
  #Print on file
  clumpfile.write(header_info)
  hfmt='#%3s '+6*'%8s '+'%10s '+'\n'
  clumpfile.write(hfmt % ('ID','phi_p','theta_p','phi_c','theta_c','dphi','dtheta','peak_cts'))
  #-------------------------------------------------------------------------------------------------------------
  #Print on screen
  print header_info
  print '# NC=%d clumps saved ' % (xpix.size)
  for kk in range(pid.size): print '# Clump oldID=%3d, newID=%3d, Ncts=%d' % (pid[kk],newid[kk],cheight[kk])
  print '#----------------------------------------------------------------'


  #Convert pix coords to physical units 
  pix_convert=xypix_converter(m,npix=npix,rangex=(xo,xf),rangey=(yo,yf))
  #Convert peak coords
  xpeak,ypeak,phipeak,thetapeak=pix_convert.get_phys_from_pix(xpix,ypix)
  #Convert peak centroid coords
  xpeakc,ypeakc,phipeakc,thetapeakc=pix_convert.get_phys_from_pix(xpix_c,ypix_c)
  #Convert peak sizes
  phi_plus_dphi,theta_minus_dtheta=pix_convert.get_phys_from_pix(xpix+dx_pix,ypix-dy_pix)[2:]
  dphi=phi_plus_dphi-phipeak
  dtheta=thetapeak-theta_minus_dtheta

  #Print peak data on file---------------------------------
  fmt='%4d '+6*'%8.3f '+'%10.0f '
  scipy.savetxt(clumpfile,np.array([newid,phipeak,thetapeak,phipeakc,thetapeakc,dphi,dtheta,cheight]).T,fmt=fmt)

  #Plot detected clump peaks
  if args.labels:
    #Peak ID labels
    m.scatter(xpeak,ypeak,c='w',alpha=0.5,edgecolor='k',s=110,zorder=100)
    for ii in range(newid.size): ax.text(xpeak[ii],ypeak[ii],newid[ii],fontsize=7,color='black',
                                        horizontalalignment='center',verticalalignment='center',zorder=101)
  else:
    m.scatter(xpeak,ypeak,c='none',edgecolor='k',s=20,zorder=99)

  #Save current figure 
  fig.savefig(figname)

  #If flag is set, plot new figure indicating pixels associated to each clump
  if args.saveclumps:
    ##Convert clump pixel mask to fits
    #print 'Converting pixel map from NDF to fits'
    #os.system('%s/convert/ndf2fits _zp_cmask _zp_cmask.fits' % (starlink_path)) 

    ##Read in fits image
    #cmask_dat=pyfits.getdata('_zp_cmask.fits',0)
    #inds=np.arange(npix)+1
    #xpix_2d,ypix_2d=np.meshgrid(inds,inds)
    #xinds,yinds=xpix_2d.flatten(),ypix_2d.flatten()
    #pcts_1d=zi.flatten()
    #cmask_1d=cmask_dat.flatten()

    #xcmask,ycmask,phicmask,thetacmask = pix_convert.get_phys_from_pix(xinds,yinds)
 
    #Plot identified clumps on top of pole count map and print out
    cmapp=plt.cm.gist_ncar(np.linspace(0., 0.9, pid.size))  #Upper limit is 0.85 to avoid last colors of the colormap
    cmapp=plt.cm.gist_ncar_r(np.linspace(0.1, 0.9, pid.size))  #Upper limit is 0.85 to avoid last colors of the colormap
    if pid.size<=10:
     cmapp=['darkviolet','orange','lime','royalblue','orchid','red','gray','pink','limegreen','navy']
#     cmapp=['orchid','red','mediumblue','orange','red','royalblue','gray','pink','limegreen','navy']

    file_clumppixfname=open(clumppixfname,'w')
    file_clumppixfname.write('#%6s %10s %10s\n' % ('IDpole','phi_pole','theta_pole'))
    for kk in np.arange(pid.size):
      #Save only pixels inside the FWXM of the peak and with counts>minheight
      pmask = (cmask_1d==pid[kk]) & (pcts_1d>=args.fwxm*cheight[kk]) & (pcts_1d>=minheight)
      #plot current peak only
      m.plot(xcmask[pmask],ycmask[pmask],color=cmapp[kk],mec='None',ms=5,marker='o',alpha=args.alpha)
      newid_cmask=newid[kk]*np.ones_like(cmask_1d[pmask])
      scipy.savetxt(file_clumppixfname,np.array([newid_cmask,phicmask[pmask],thetacmask[pmask]]).T,fmt='%7d %10.4f %10.4f')

    cfigname='%s.%s.%s.%s.pls.%s' % (figname_root,mode,proj[:3],pmode,args.fig)
    fig.savefig(cfigname)

    #Remove auxiliary files
    os.system('rm -f _zp* _zp.fits')

  if args.show: plt.show()
  else: fig.clf()
