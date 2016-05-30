#!/usr/bin/env python
import pylab as plt
import scipy
import scipy.interpolate
import numpy as np
from mpl_toolkits.basemap import Basemap
import sys
import argparse
import scipy.ndimage
import scipy.interpolate
import pyfits
import os
import myutils

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
parser.add_argument('-m',help='Plot mGC3/nGC3/GC3/mGC3hel/AUX pole count map. Default is mGC3', action='store',default='mGC3',choices=['mGC3','nGC3','GC3','AUX','mGC3hel'])
parser.add_argument('-f','--fig',help='Output plot type png/eps. Default is png', action='store',default='png',choices=['png','eps','pdf'])
parser.add_argument('-ext',metavar='outfile_ext',help='Output suffix [optional]. If given output will be infile.outfile_ext.mgc3.pst',action='store',default=['',],nargs=1)
parser.add_argument('-log',help='Plot detected peaks in log-count map', action='store_true',default=False)
parser.add_argument('-nolabels',help='Plot peak ID labels', action='store_true',default=False)
parser.add_argument('-title',help='Plot title', action='store',default=None)
parser.add_argument('-lon0',help='Longitude for Y-axis. Default is 0.', action='store',default=0.,type=np.float)
parser.add_argument('-lat0',help='Bounding latitude for plot. Default is 90.', action='store',default=0.,type=np.float)
parser.add_argument('-latmax',help='Max latitude upto which meridians are drawn. Default is 80.', action='store',default=80.,type=np.float)
parser.add_argument('-mlab','--merlabels',help='Show meridian labels. Default is False', action='store_true',default=False)
parser.add_argument('-mlabr','--merlabelsr',help='Show meridian labels (right axes). Default is False', action='store_true',default=False)
parser.add_argument('-plab','--parlabels',help='Show parallel labels. Default is False', action='store_true',default=False)
parser.add_argument('-vmin',help='Min counts for color-scale. Default is min(cts)', action='store',default=None,type=np.float)
parser.add_argument('-vmax',help='Max counts for color-scale. Default is max(cts)', action='store',default=None,type=np.float)
parser.add_argument('-dlat',help='Spacing between parallels. Default is 20.', action='store',default=20.,type=np.float)
parser.add_argument('-dlon',help='Spacing between meridians. Default is 30.', action='store',default=30.,type=np.float)
parser.add_argument('-ms',help='Marker size. Default: 50 for npaeqd.', action='store',default=50,type=np.float)
parser.add_argument('-r','--raw',help='Plot raw grid pole-count map instead of contour map.', action='store_true',default=False)
parser.add_argument('-t','--twohemispheres',help='Plot both hemispheres in pole-count map.', action='store_true',default=False)
parser.add_argument('-force_onehemisph',help='Force using one hemisphere allways', action='store_true',default=False)
parser.add_argument('-s','--show',help='Show plot in window. Default is False', action='store_true',default=False)
parser.add_argument('-nc','--noclumps',help='Do not plot or save poles associated to each peak.', action='store_true',default=False)
parser.add_argument('-bw',help='Use grayscale colormap to plot PCMs. Default False (uses jet colormap)', action='store_true',default=False)
parser.add_argument('-cmap',help='Choose color map. Default is sron', action='store',default='sron',choices=['sron','gray','gray_r','viridis','inferno'])
parser.add_argument('-npix',help='Default number of pixels to resample image for contour plotting', action='store',type=np.int,default=500)
parser.add_argument('-mj','--maxjump',help='Fellwalker MaxJump param, neighbourhood radius to search for +gradient. Default 6.', action='store',default=6,type=np.float)
parser.add_argument('-md','--mindip',help='Fellwalker MinDip param, two clumps are merged if height difference <MinDip. Default 2*RMS', action='store',default=None)
parser.add_argument('-al','--alpha',help='Clump transparency. Default 0.4', action='store',default=0.4,type=np.float)
peakargs = parser.add_mutually_exclusive_group()
peakargs.add_argument('-fr','--frms',help='Default option. Min peak height is frms*RMS. Default fr=5.', action='store',type=np.float,default=5.)
peakargs.add_argument('-ff','--ffrac',help='Min peak height is fmax*max_pole_counts', action='store',type=np.float)
parser.add_argument('-fx','--fwxm',help='Store pixels with cts>fwxm*maxpeak. Default is 0.5. If -U, fx is in N-sigma units', action='store',default=0.5,type=np.float)
parser.add_argument('-flab','--flabels',help='Increase size of peak labels by factor flab. Default 1.', action='store',default=1.0,type=np.float)
parser.add_argument('-fcirc','--fcirc',help='Increase size of peak markers by factor fcirc. Default 1.', action='store',default=1.0,type=np.float)
parser.add_argument('-U','--unsharp',help='Detect peaks in unsharp-masked image. Default is False', action='store_true',default=False)
parser.add_argument('-ns','--nsigma',help='If -U is set, set N-sigma threshold for detections. Default is 3.', action='store',type=np.float,default=3.)
parser.add_argument('-nmed',help='If -U is set, size of neighbourhood for median computation. Default is 60.', action='store',type=np.float,default=60.)


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
if 'hel' in mode:   counts_col=4-1
elif 'mgc3' in mode:   counts_col=3-1
elif 'ngc3' in mode: counts_col=6-1
elif 'gc3'  in mode: counts_col=5-1
elif 'aux'  in mode: counts_col=8-1

#Parse raw/contour mode-------------------------
if args.raw: 
  pmode='r'
  print 'Plotting raw pole-count map'
else: 
  pmode='c'
  print 'Plotting contour pole-count map'
#If log-flag is set, do everything with log(counts)
if args.log: 
   print 'Do peak detection on linear scale but show plot on log scale'
   pmode=pmode+'l'

if args.twohemispheres:
  print 'Plotting both hemispheres in pole-count map'
else: print 'Plotting one hemisphere in pole-count map.'

proj='npaeqd'
print 'Plotting using projection:', proj

ori='vertical'
ori='horizontal'
ni=0

#colormap=plt.cm.jet
if 'inferno' in args.cmap:
 colormap=newcmap.inferno
elif 'viridis' in args.cmap:
 colormap=newcmap.viridis
elif 'gray' in args.cmap:
 if '_r' not in args.cmap: colormap=plt.cm.gray
 else: colormap=plt.cm.gray_r
else:
 colormap=myutils.get_sron_rainbow(N=11)

colormap_nsig=plt.cm.spectral

for infilen in file_list:

  phio,thetao,pole_ctso=pdat=scipy.genfromtxt(infilen,comments='#',usecols=(0,1,counts_col),unpack=True)

  #Default title----------------------------------
  if not args.title: args.title=infilen

  #Output figure and file names
  figname_root=infilen.replace('.mgc3.cts',args.ext[0])  #works well if args.ext is empty
  figname='%s.%s.%s.%s.%s' % (figname_root,mode,proj[:3],pmode,args.fig)
  clumpfname='%s.%s.peak.pls' % (figname_root,mode)
  clumppixfname='%s.%s.pls' % (figname_root,mode)
  print 'Output filename:', figname

  phi2=np.append(phio,(phio+180.) % 360.)
  theta2=np.append(thetao,-thetao)
  pole_cts2=np.append(pole_ctso,pole_ctso) 

  if args.twohemispheres or args.unsharp and not args.force_onehemisph:
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
  grid_color=(0.65,0.65,0.65)
  mlabels_dic={'labels':[0,0,0,0],'labelstyle':'+/-'}
  if args.merlabels:  mlabels_dic['labels'][0]=1
  if args.merlabelsr: mlabels_dic['labels'][1]=1
  if args.parlabels: plabels_dic={'labels':[0,0,0,1],'labelstyle':'+/-'}
  else: plabels_dic={'labels':[0,0,0,0]}
  m.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]),color=grid_color,linewidth=1.,
                   latmax=args.latmax,**mlabels_dic)
  m.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]),color=grid_color,linewidth=1.,**plabels_dic)
  m.drawmapboundary()
  x,y=m(phi,theta)

  #------------Grid-data for contour plotting---------------
  #npix=500
  npix=args.npix
  clevels=30
  xo,xf=np.min(x),np.max(x)
  yo,yf=np.min(x),np.max(x)
  xi = np.linspace(xo,xf,npix)
  yi = np.linspace(yo,yf,npix)
  #1-D flattened arrays
  xi1d,yi1d=np.meshgrid(xi,yi)
  pi1d,ti1d=m(xi1d,yi1d,inverse=True)
  zi1 = scipy.interpolate.griddata((x,y), pole_cts, (xi1d,yi1d), method='linear')
  tmask=(ti1d<-90.) | (ti1d>90.)
  zi=np.ma.masked_array(zi1,mask=tmask) #These masks in masked-arrays work the other way, True are masked-OUT values

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
   zii=zi.copy()
   if args.log: 
     if not args.twohemispheres: zii[ti1d<0.]=0.
     c=m.contourf(xi,yi,np.log10(zii),clevels,cmap=colormap,vmin=args.vmin,vmax=args.vmax,antialiased=False)
   else:
     if not args.twohemispheres: zii[ti1d<0.]=np.nan
     c=m.contourf(xi,yi,zii/10**lmax, clevels,cmap=colormap,vmin=args.vmin,vmax=args.vmax)

  #Plot colorbar
  cax0=plt.gca().get_position()
  cax=plt.axes([cax0.x0,cax0.y0+dw+0.05,dw,0.02])
  cb=plt.colorbar(c,cax=cax,orientation='horizontal',format='%4.1f')
  cax.xaxis.set_ticks_position('top')
  
  #Labels and such
  if lmax>0: factorl='$\\times 10^{%d}$ ' % (lmax)
  else: factorl=''
  if args.log:
     cax.set_xlabel('%s log-pole-counts (dex stars/pole)' % (mode_ori))
  else:
     cax.set_xlabel('%s pole-counts (%sstars/pole)' % (mode_ori,factorl))
  cax.xaxis.set_label_position('top') 

  if args.title:
    ax.text(0.5,1.14,args.title,transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontsize=16)

  #-------------------------Unsharp masking-----------------------------------------------------
  if args.unsharp:
    #Compute median filtered image
    nsm=args.nmed #neighbourhood for median computation
    zi_smooth=scipy.ndimage.median_filter(zi,size=(nsm,nsm),mode='wrap')
    zi_sharp=zi-zi_smooth
    zi_sigma=np.sqrt(zi_smooth)  #assuming counts in the PCM follow a Poisson distribution
    zi_sharp_Nsig=zi_sharp/zi_sigma   #express subtracted image in nsigma-units
    #----plot unsharp masked and smoothed images
    fig2=plt.figure(2,figsize=(16,6.3))
    fig2.subplots_adjust(bottom=0.01,top=0.99,left=0.03,right=0.97,wspace=0.1)
    #----------------smoothed image---------------------------------
    axs=fig2.add_subplot(1,3,1)
    ms = Basemap(projection=proj,ax=axs,**proj_dict)
    ms.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]),color='grey',lw=2.)
    ms.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]),color='grey',lw=2.)
    ms.drawmapboundary()
    if not args.twohemispheres: zi_smooth[ti1d<0.]=np.nan  #Unless -t is explicitly set, don't plot S-counts
    c1=ms.contourf(xi,yi,zi_smooth, clevels,cmap=colormap,vmin=0.)
    cb=plt.colorbar(c1,ax=axs,orientation='horizontal',format='%d',pad=0,aspect=30)
    cb.set_label('%s pole-counts (stars/pole)' % (mode_ori),fontsize=15)
    axs.set_title('Smoothed %s PCM' % (mode_ori),fontsize=15)
    #-------------subtracted image-----------------------------------------
    axu=fig2.add_subplot(1,3,2)
    mu = Basemap(projection=proj,ax=axu,**proj_dict)
    mu.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]),color='lightgrey',lw=2.)
    mu.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]),color='lightgrey',lw=2.)
    mu.drawmapboundary()
    if not args.twohemispheres: zi_sharp[ti1d<0.]=0.  #Unless -t is explicitly set, don't plot S-counts
    c3=mu.contourf(xi,yi,np.log10(zi_sharp),clevels,cmap=colormap,vmin=0.)
    cb=plt.colorbar(c3,ax=axu,orientation='horizontal',format='%4.1f',pad=0,aspect=30)
    if args.log:
       cb.set_label('%s log-pole-counts (dex stars/pole)' % (mode_ori),fontsize=15)
    else:
       cb.set_label('%s pole-counts (%sstars/pole)' % (mode_ori,factorl),fontsize=15)
    axu.set_title('Unsharp-masked %s PCM' % (mode_ori),fontsize=15)
    #----------------subtracted image in Nsigma units---------------------------------
    axn=fig2.add_subplot(1,3,3)
    mn = Basemap(projection=proj,ax=axn,**proj_dict)
    mn.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]),color='lightgrey',lw=2.)
    mn.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]),color='lightgrey',lw=2.)
    mn.drawmapboundary()
    zi_sharp_cut=zi_sharp_Nsig.copy()
    sigmax=12.  #Maximum Nsigma for contour and color display
    zi_sharp_cut[zi_sharp_cut>=sigmax]=sigmax
    if not args.twohemispheres: zi_sharp_cut[ti1d<0.]=np.nan  #Unless -t is explicitly set, don't plot S-counts
    c2=mn.contourf(xi,yi,(zi_sharp_cut), np.arange(0.,sigmax+1.,1.),cmap=colormap_nsig)
    cb=plt.colorbar(c2,ax=axn,orientation='horizontal',format='%4.1f',pad=0,aspect=30,extend='max')
    #Labels and such
    cb.set_label('Significance ($N\sigma$)',fontsize=15)
    axn.set_title('Unsharp-masked %s PCM (N$\sigma$)' % (mode_ori),fontsize=15)
    #---------
    axu.text(0.5,0.97,args.title,transform=fig2.transFigure,horizontalalignment='center',
             verticalalignment='center',fontsize=17)
    usharp_figname='%s.%s.%s.%s.usharp.%s' % (figname_root,mode,proj[:3],pmode,args.fig)
    #fig2.set_rasterized(True)
    fig2.savefig(usharp_figname)
    #-------------------------------Print smoothed and unsharp masked PCMs------------
    usharp_filen='%s.%s.%s.%s.usharp.cts' % (figname_root,mode,proj[:3],pmode)
    head='%3s %7s %8s %8s %8s' % ('phi','theta','Nsmooth','Nusharp','Nusharp_Nsig')
    pi1d,ti1d,zi_smooth_1d=pi1d.flatten(),ti1d.flatten(),zi_smooth.flatten()
    pi1d[pi1d<0]=pi1d[pi1d<0]+360.
    mask=(ti1d>0) & (ti1d<=90) & (zi_smooth.flatten()>0)
    usharp_file=open(usharp_filen,'w')
    usharp_file.write('#Smoothed map stats: \n')
    usharp_file.write('# Min=%.0f Max=%.0f Mean=%.0f StdDev=%.0f\n# Median=%.0f P10=%.0f P90=%.0f\n' 
                      % (np.min(zi_smooth_1d[mask]),np.max(zi_smooth_1d[mask]),
                         np.mean(zi_smooth_1d[mask]),np.std(zi_smooth_1d[mask]),
                         np.median(zi_smooth_1d[mask]),np.percentile(zi_smooth_1d[mask],10),
                         np.percentile(zi_smooth_1d[mask],90)))
    scipy.savetxt(usharp_file,np.array([pi1d[mask],ti1d[mask],zi_smooth.flatten()[mask],zi_sharp.flatten()[mask],
                  zi_sharp_Nsig.flatten()[mask]]).T,fmt='%6.3f %7.3f %8.1f %8.1f %8.1f',
                  header=head)
    #pi1d,ti1d,zi_smooth_1d=pi1d[mask],ti1d[mask],zi_smooth_1d[mask]

  if args.noclumps:
    #fig.set_rasterized(True)
    fig.savefig(figname)
    if args.unsharp: fig2.savefig(usharp_figname)
    if args.show: plt.show()
    else: fig.clf()
    print 'WARNING: No-clump flag set, plain PCM plotted. Skipping peak detection...  '
    continue


  #--------------------Manage formats to run the Fellwaker peak detection algorithm---------------------------------------------
  if args.unsharp: zi=zi_sharp_Nsig
 
  hdu = pyfits.PrimaryHDU()
  print np.shape(zi), np.shape(zi.data)
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
     rms=1.  #If -U, zi is always in N-sigma units by default            
     minheight=args.nsigma
  else:
   rms=np.std(zi)
   if args.frms is not None:
     minheight=args.frms*rms
   else: 
     minheight=args.ffrac*np.max(zi)
  if args.mindip is None:
    mindip='' #Will use Fellwalker's default
  else:
    mindip='config="fellwalker.mindip=%f" ' % (np.float(args.mindip))
  print 'Finding clumps with Fellwalker'
  os.system('%s/cupid/findclumps in=_zpndf.sdf out=_zp_cmask method=fellwalker outcat=_zp_clumps rms=%f config="fellwalker.maxjump=%.0f" %s ' % (starlink_path,rms,args.maxjump,mindip)) 

  #------------------------------Deal with identified peaks (peak centroids, etc)---------------------------
  #Read-in output table with identified clumps
  clumpdat=pyfits.open('_zp_clumps.FIT')[1].data                    
  xpix,ypix,cheight=clumpdat.field('Peak1'),clumpdat.field('Peak2'),clumpdat.field('Peak')
  xpix_c,ypix_c,dx_pix,dy_pix=clumpdat.field('Cen1'),clumpdat.field('Cen2'),clumpdat.field('Size1'),clumpdat.field('Size2')
  pid=np.arange(xpix.size) + 1  #Peak IDs, must start from 1

  #Keep only peaks with counts > minheight and with sizes>=1 (clumps have to be at least 1 pixel in size)
  mask=(cheight>=minheight) #& ((dx_pix>=1.) | (dy_pix>=1.)) 
  #mask=(cheight>=minheight) & ((dx_pix>=1.) | (dy_pix>=1.)) 
  #mask=(cheight>=minheight) & ((dx_pix>=2.) | (dy_pix>=2.)) 
  if not mask.any(): 
    print 'No peaks above minheight or with size>=1pix found'
    continue
  xpix,ypix,cheight=xpix[mask],ypix[mask],cheight[mask]
  xpix_c,ypix_c,dx_pix,dy_pix,pid=xpix_c[mask],ypix_c[mask],dx_pix[mask],dy_pix[mask],pid[mask]

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
  if args.unsharp: smooth_cts_1d,sharp_cts_1d=zi_smooth_1d,zi_sharp.flatten()
  cmask_1d=cmask_dat.flatten() #cmask goes with pid, i.e. pid=cmask_1d.unique()
  xcmask,ycmask,phicmask,thetacmask = pix_convert.get_phys_from_pix(xinds,yinds)
  #----------Deal with clumps that go from one hemisphere to the other-------------------------------
  #Loop over clump IDs, if centroid is in the N, save all its pixels
  #Peak data
  u_phipeak,u_thetapeak,u_xpix,u_ypix=np.array([]),np.array([]),np.array([]),np.array([])
  u_dx_pix,u_dy_pix,u_pid,u_cheight=np.array([]),np.array([]),np.array([]),np.array([])
  u_xpix_c,u_ypix_c=np.array([]),np.array([])
  #Pixel data
  u_xcmask,u_ycmask=np.array([]),np.array([])
  u_phicmask,u_thetacmask=np.array([]),np.array([])
  u_cmask_1d,u_pcts_1d=np.array([]),np.array([])
  if args.unsharp: u_smooth_cts_1d,u_sharp_cts_1d=np.array([]),np.array([])
  for ii in range(thetapeak.size):
   if thetapeak[ii]>=0: #North peaks
     #Save peak data
     u_phipeak,u_thetapeak=np.append(u_phipeak,phipeak[ii]),np.append(u_thetapeak,thetapeak[ii])
     u_xpix,u_ypix=np.append(u_xpix,xpix[ii]),np.append(u_ypix,ypix[ii])   
     u_xpix_c,u_ypix_c=np.append(u_xpix_c,xpix_c[ii]),np.append(u_ypix_c,ypix_c[ii])   
     u_dx_pix,u_dy_pix=np.append(u_dx_pix,dx_pix[ii]),np.append(u_dy_pix,dy_pix[ii])   
     u_pid,u_cheight=np.append(u_pid,pid[ii]),np.append(u_cheight,cheight[ii])
     #Save pixel-peak data
     peakmask=(cmask_1d==pid[ii]) 
     u_xcmask,u_ycmask=np.append(u_xcmask,xcmask[peakmask]),np.append(u_ycmask,ycmask[peakmask])
     u_phicmask,u_thetacmask=np.append(u_phicmask,phicmask[peakmask]),np.append(u_thetacmask,thetacmask[peakmask])
     u_cmask_1d,u_pcts_1d=np.append(u_cmask_1d,cmask_1d[peakmask]),np.append(u_pcts_1d,pcts_1d[peakmask])
     if args.unsharp: 
      u_smooth_cts_1d=np.append(u_smooth_cts_1d,smooth_cts_1d[peakmask])
      u_sharp_cts_1d =np.append(u_sharp_cts_1d,sharp_cts_1d[peakmask])
   
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
        mask_tol=dist2d.deg<0.5 #Keep only matches within less than 1deg
        cmaskc[mask_tol]=u_cmask_1d[Nindex[mask_tol]]
        #Put everything back together 
        if mask_tol.any():
          mask_tol=(dist2d.deg<0.5) & (thetac<0) #Keep only matches within less than 1deg
          #Save Pixel data (these are, by definition, repeated peaks so nothing should be appended to peak-data arrays)
          u_xcmask,u_ycmask=np.append(u_xcmask,xcmask[peakmask][mask_tol]),np.append(u_ycmask,ycmask[peakmask][mask_tol])
          u_phicmask,u_thetacmask=np.append(u_phicmask,phicmask[peakmask][mask_tol]),np.append(u_thetacmask,thetacmask[peakmask][mask_tol])
          u_cmask_1d,u_pcts_1d=np.append(u_cmask_1d,cmaskc[mask_tol]),np.append(u_pcts_1d,pcts_1d[peakmask][mask_tol])
          if args.unsharp: 
           u_smooth_cts_1d=np.append(u_smooth_cts_1d,smooth_cts_1d[peakmask][mask_tol])
           u_sharp_cts_1d =np.append(u_sharp_cts_1d,sharp_cts_1d[peakmask][mask_tol])
     else:
       continue

  #Dump any pixels leftover in the S
  tmask=u_thetacmask>=0.
  u_xcmask,u_ycmask=u_xcmask[tmask],u_ycmask[tmask]
  u_phicmask,u_thetacmask=u_phicmask[tmask],u_thetacmask[tmask]
  u_cmask_1d,u_pcts_1d=u_cmask_1d[tmask],u_pcts_1d[tmask]
  if args.unsharp: 
   u_smooth_cts_1d,u_sharp_cts_1d=u_smooth_cts_1d[tmask],u_sharp_cts_1d[tmask]

  #u_newid=np.arange(u_phipeak.size) + 1  #Rename so IDs will be consecutive numbers starting from 1

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
    if args.mindip is not None: header_info=header_info+'#         MinDip=%-4.1f (Nsigma-units)\n' % (np.float(args.mindip))
  else:
   if args.frms is not None: 
     header_info=header_info+'#         MinHeight=%d (=%.1f*RMS)\n' % (minheight,args.frms)
   else: 
    header_info=header_info+'#         MinHeight=%d (=%.2f*max_cts)\n' % (minheight,args.ffrac)
   if args.mindip is not None: header_info=header_info+'#         MinDip=%-4.1f (cts)\n' % (np.float(args.mindip))
  if args.unsharp: header_info=header_info+'#         Nmedian=%d\n' % (args.nmed)
  header_info=header_info+'#         FWXM=%.4f (stored pixels associated to each clump)\n' % (args.fwxm)
  header_info=header_info+'#----------------------------------------------------------------------\n'
  #Print on file
  clumpfile.write(header_info)
  hfmt='#%3s '+6*'%8s '+'%10s %10s %10s %10s'+'\n'
  clumpfile.write(hfmt % ('ID','phi_p','theta_p','phi_c','theta_c','dphi','dtheta','peak_cts','exp_purity','Nsmooth','Nsharp'))
  #-------------------------------------------------------------------------------------------------------------
  #Print on screen
  print header_info
  print '# NC=%d clumps saved ' % (u_xpix.size)
  print '#----------------------------------------------------------------'

  ##Convert pix coords to physical units 
  ##Convert peak coords
  u_xpeak,u_ypeak,u_phipeak,u_thetapeak=pix_convert.get_phys_from_pix(u_xpix,u_ypix)
  ##Convert peak centroid coords
  u_xpeakc,u_ypeakc,u_phipeakc,u_thetapeakc=pix_convert.get_phys_from_pix(u_xpix_c,u_ypix_c)
  #Convert peak sizes
  phi_plus_dphi,theta_minus_dtheta=pix_convert.get_phys_from_pix(u_xpix+u_dx_pix,u_ypix-u_dy_pix)[2:]
  dphi=phi_plus_dphi-u_phipeak
  dtheta=u_thetapeak-theta_minus_dtheta

  #Compute expected purity
  #exp_purity,Nsm_l,Nsh_l=np.array([]),np.array([]),np.array([])
  fmt='%4d '+6*'%8.3f '+'%10.1f %10.3f %10.0f %10.0f\n'
  newid,u_newid=0,np.array([])

  for kk in range(u_pid.size):  
    #Using only the pixels that will be printed out
    imask = (u_cmask_1d==u_pid[kk]) 
    if not imask.any(): continue
    if args.unsharp:
     pmask = (u_cmask_1d==u_pid[kk]) & (u_sharp_cts_1d>=args.fwxm*np.max(u_sharp_cts_1d[imask])) & (u_pcts_1d>=minheight)
     if pmask.any():
      Nsm,Nsh=np.sum(u_smooth_cts_1d[pmask]),np.sum(u_sharp_cts_1d[pmask])
      Nsigma=np.max(u_sharp_cts_1d[pmask].astype(float)/np.sqrt(u_smooth_cts_1d[pmask]))
      exp_purity=Nsh/np.float(Nsh+Nsm)
     else:
      Nsm,Nsh=0.,0.
      Nsigma,exp_purity=0.,-1
    else:
     pmask = (u_cmask_1d==u_pid[kk]) & (u_pcts_1d>=args.fwxm*u_cheight[kk]) & (u_pcts_1d>=minheight)
     Nsm,Nsh=0,0
     exp_purity=-1
    #if pmask.any(): 
    if pmask.sum()>5: #Require at least five pixels above threshold
     newid=newid+1
     u_newid=np.append(u_newid,newid)
     #Fix peak height, this was wrong because Fellwalker returns the sum over each clump (I think)
     if args.unsharp: u_cheight[kk]=Nsigma
     print '# Clump oldID=%3d, newID=%3d, height=%.1f' % (u_pid[kk],newid,u_cheight[kk])
     #Print peak data on file---------------------------------
     clumpfile.write(fmt % (newid,u_phipeak[kk],u_thetapeak[kk],u_phipeakc[kk],u_thetapeakc[kk],
                            dphi[kk],dtheta[kk],u_cheight[kk],exp_purity,Nsm,Nsh))
    else: 
     u_newid=np.append(u_newid,0)
     print '# Skipping clump oldID=%3d, less than 5 pixels>threshold_height'  % (u_pid[kk])
  print '#----------------------------------------------------------------'

  if args.noclumps:  
    #fig.set_rasterized(True)
    fig.savefig(figname)

  #Plot detected clump peaks and labels
  if args.unsharp: axl,ml=[ax,axn,axu],[m,mn,mu]
  else: axl,ml=[ax,],[m,]
  if not args.nolabels:
    #Peak ID labels
    for aax,mm in zip(axl,ml):
      mm.scatter(u_xpeak[u_newid>0],u_ypeak[u_newid>0],c='w',alpha=0.5,edgecolor='k',s=110*args.fcirc,zorder=100)
      for ii in range(u_newid.size):
       if u_newid[ii]>0: 
         aax.text(u_xpeak[ii],u_ypeak[ii],'%d' % (u_newid[ii]),fontsize=7*args.flabels,color='black',
                   horizontalalignment='center',verticalalignment='center',zorder=101)
  else:
    for mm in ml:
      mm.scatter(u_xpeak,u_ypeak,c='none',edgecolor='k',s=20,zorder=99)

  #Save current figure after plotting peaks and labels
  if not args.noclumps: 
    fig.savefig(figname)
    if args.unsharp: fig2.savefig(usharp_figname)

  #If flag is set, plot new figure indicating pixels associated to each clump
  if not args.noclumps:
    #Plot identified clumps on top of pole count map and print out
    #cmapp=plt.cm.gist_ncar(np.linspace(0., 0.9, u_pid.size))  #Upper limit is 0.85 to avoid last colors of the colormap
    cmapp=plt.cm.gist_ncar_r(np.linspace(0.1, 0.9, u_pid.size))  #Upper limit is 0.85 to avoid last colors of the colormap
    if u_pid.size<=10:
     cmapp=['darkviolet','orange','lime','royalblue','orchid','red','gray','pink','limegreen','navy']
#     cmapp=['orchid','red','mediumblue','orange','red','royalblue','gray','pink','limegreen','navy']

    file_clumppixfname=open(clumppixfname,'w')
    file_clumppixfname.write('#%6s %10s %10s %10s %10s\n' % ('IDpole','phi_pole','theta_pole','Nsmooth','Nsharp'))

    for kk in np.arange(u_newid.size):
     if u_newid[kk]>0:
      #Save only pixels inside the FWXM of the peak and with counts>minheight
      imask = (u_cmask_1d==u_pid[kk]) 
      if args.unsharp:
       pmask = (u_cmask_1d==u_pid[kk]) & (u_sharp_cts_1d>=args.fwxm*np.max(u_sharp_cts_1d[imask])) & (u_pcts_1d>=minheight)
      else:
       pmask = (u_cmask_1d==u_pid[kk]) & (u_pcts_1d>=args.fwxm*u_cheight[kk]) & (u_pcts_1d>=minheight)
      #plot current peak only
      m.scatter(u_xcmask[pmask],u_ycmask[pmask],c=cmapp[kk],edgecolors='none',s=20,marker='o',alpha=args.alpha)
      u_newid_cmask=u_newid[kk]*np.ones_like(u_cmask_1d[pmask])
      #Sum smoothed counts over each clump
      if args.unsharp: Nsm,Nsh=u_smooth_cts_1d[pmask],u_sharp_cts_1d[pmask]
      else: Nsm,Nsh=np.zeros_like(u_phicmask[pmask]),np.zeros_like(u_phicmask[pmask])
      scipy.savetxt(file_clumppixfname,np.array([u_newid_cmask,u_phicmask[pmask],u_thetacmask[pmask],Nsm,Nsh]).T,
                    fmt='%7d %10.4f %10.4f %10.1f %10.1f')

    cfigname='%s.%s.%s.%s.pls.%s' % (figname_root,mode,proj[:3],pmode,args.fig)
    fig.savefig(cfigname)

    #Remove auxiliary files
    os.system('rm -f _zp* _zp.fits')
  
  if args.show: plt.show()
  else: fig.clf()
  fig.clf()
  if args.unsharp: fig2.clf()
