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
    phic=phic+360.

    return (xc,yc,phic,thetac)

#----------------------------------------------
parser = argparse.ArgumentParser(description='Detect peaks in mGC3/nGC3/GC3 pole count maps')
parser.add_argument('infile',metavar='infile',help='Input file containing pole count maps (*.cts file)',nargs=1,action='store')
parser.add_argument("-l", "--llist", action="store_true",help='Take infile as list of mgc3.cts files')
parser.add_argument('-m',help='Plot mGC3/nGC3/GC3 pole count map. Default is mGC3', action='store',default='mGC3',choices=['mGC3','nGC3','GC3'])
parser.add_argument('-f','--fig',help='Output plot type png/eps. Default is png', action='store',default='png',choices=['png','eps','pdf'])
parser.add_argument('-lon0',help='Longitude for Y-axis. Default is 0.', action='store',default=0.,type=np.float)
parser.add_argument('-lat0',help='Bounding latitude for plot. Default is 90.', action='store',default=0.,type=np.float)
parser.add_argument('-dlat',help='Spacing between parallels. Default is 30.', action='store',default=20.,type=np.float)
parser.add_argument('-dlon',help='Spacing between meridians. Default is 30.', action='store',default=30.,type=np.float)
parser.add_argument('-ms',help='Marker size. Default: 10 for npaeqd.', action='store',default=15,type=np.float)
parser.add_argument('-c','--contour',help='Plot pole-count contour map instead of raw grid.', action='store_true',default=False)
parser.add_argument('-t','--twohemispheres',help='Plot both hemispheres in pole-count map.', action='store_true',default=False)
parser.add_argument('-s','--show',help='Show plot in window. Default is False', action='store_true',default=False)
parser.add_argument('-sc','--saveclumps',help='Plot and save poles associated to each peak.', action='store_true',default=False)
peakargs = parser.add_mutually_exclusive_group()
peakargs.add_argument('-frms',help='If set, min peak height is frms*RMS', action='store',type=np.float)
peakargs.add_argument('-ffrac',help='Default option. Min peak height is fmax*max_pole_counts. Default fmax=0.6', action='store',default=0.6,type=np.float)


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
if args.contour: 
  pmode='c'
  print 'Plotting contour pole-count map'
else: 
  pmode='r'
  print 'Plotting raw pole-count map'

if args.twohemispheres:
  print 'Plotting both hemispheres in pole-count map'
else: print 'Plotting one hemisphere in pole-count map.'

proj='npaeqd'
print 'Plotting using projection:', proj

ori='vertical'
ori='horizontal'
ni=0

colormap=plt.cm.jet
#colormap=plt.cm.spectral
for infilen in file_list:

  phio,thetao,pole_ctso=pdat=scipy.genfromtxt(infilen,comments='#',usecols=(0,1,counts_col),unpack=True)
  figname_root=infilen.replace('.mgc3.cts','')
  figname='%s.%s.%s.%s.%s' % (figname_root,mode,proj[:3],pmode,args.fig)
  clumpfname='%s.%s.pls.peak.dat' % (figname_root,mode)
  clumppixfname='%s.%s.pls' % (figname_root,mode)
  print 'Output filename:', figname

  phi2=np.append(phio,(phio+180.) % 360.)
  theta2=np.append(thetao,-thetao)
  pole_cts2=np.append(pole_ctso,pole_ctso) 

  if args.twohemispheres:
   phi,theta,pole_cts=phi2,theta2,pole_cts2
  else: 
   phi,theta,pole_cts=phio,thetao,pole_ctso 
  #This is done this way so that for contour plotting, the twohemispheres are used for the interpolation
  #into a regular grid, but only one is shown 

  #----------------Pole count map------------------------------
  mer_grid=[0.,360.,args.dlon]
  par_grid=[-90.,+90.,args.dlat]

  #For npa and moll projections, plot map as viewed from lon0 only
  fig=plt.figure(1,figsize=(8,8))
  dw=0.8
  wo=(1.-dw)/2.
  wyo=0.75*wo
  fig.subplots_adjust(left=wo,right=dw+wo,top=dw+wo,bottom=wyo)
  nrow,ncol,nplot=1,1,1
  l0=args.lon0
  proj_dict={'boundinglat':args.lat0,'resolution':'l','lon_0':l0}
  ms=args.ms

  ax=fig.add_subplot(nrow,ncol,nplot)
  m = Basemap(projection=proj,ax=ax,**proj_dict)
  m.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]),color='white',lw=2.)
  m.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]),color='white',lw=2.)
  m.drawmapboundary()

  x,y=m(phi,theta)
  #
  npix=400
  clevels=30
  xo,xf=np.min(x),np.max(x)
  yo,yf=np.min(x),np.max(x)
  xi = np.linspace(xo,xf,npix)
  yi = np.linspace(yo,yf,npix)
  zi = plt.griddata(x,y,pole_cts,xi,yi) #,'nn')

  lmax=np.floor(np.log10(np.max(pole_cts)))
  if 'r' in pmode: c=m.scatter(x,y,c=pole_cts/10**lmax,edgecolor='none',s=ms,cmap=colormap)
  else: c=m.contourf(xi,yi,zi/10**lmax,clevels,cmap=colormap)

  #Plot colorbar
  cax0=plt.gca()
  cax=plt.axes([wo,1.2*wyo+dw,dw,0.02])
  cb=plt.colorbar(c,cax=cax,orientation='horizontal',format='%4.1f',label='prueba')
  cax.xaxis.set_ticks_position('top')
  
  #Labels and such
  if lmax>0: factorl='$\\times 10^{%d}$ ' % (lmax)
  cax.set_xlabel('%s pole-counts (%sstars/pole)' % (mode_ori,factorl))
  cax.xaxis.set_label_position('top') 

  #print fits image
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
  rms=np.sqrt(np.mean(zi))
  if args.frms is not None:
    minheight=args.frms*rms
  else: 
    minheight=args.ffrac*np.max(zi)
  print 'Finding clumps with Fellwalker'
  os.system('%s/cupid/findclumps in=_zpndf.sdf out=_zp_cmask method=fellwalker outcat=_zp_clumps rms=%f' 
            % (starlink_path,rms)) 

  #Read-in output table with identified clumps
  clumpdat=pyfits.open('_zp_clumps.FIT')[1].data                    
  xpix,ypix,cheight=clumpdat.field('Peak1'),clumpdat.field('Peak2'),clumpdat.field('Peak')
  xpix_c,ypix_c,dx_pix,dy_pix=clumpdat.field('Cen1'),clumpdat.field('Cen2'),clumpdat.field('Size1'),clumpdat.field('Size2')
  pid=np.arange(xpix.size) + 1  #Peak IDs, must start from 1

  #Keep only peaks with counts > minheight
  mask=cheight>=minheight
  if not mask.any(): continue
  xpix,ypix,cheight=xpix[mask],ypix[mask],cheight[mask]
  xpix_c,ypix_c,dx_pix,dy_pix,pid=xpix_c[mask],ypix_c[mask],dx_pix[mask],dy_pix[mask],pid[mask]

  #print pars on screen
  print '#----------------------------------------------------------------'
  print '# Peak-detection algorithm: Starlink Fellwalker (Berry+2014)'
  print '# Params: RMS=%.1f (=sqrt(mean(cts))' % (rms)
  if args.frms is not None:
    print '#         MinHeight=%d (=%.1f*RMS)' % (minheight,args.frms)
  else:
    print '#         MinHeight=%d (=%.2f*max_cts)' % (minheight,args.ffrac)
  print '#----------------------------------------------------------------'
  print '# NC=%d clumps saved ' % (xpix.size)
  for kk in range(pid.size): print '# Clump N=%d, Ncts=%d' % (pid[kk],cheight[kk])
  print '#----------------------------------------------------------------'

  #Open output file to store clump coords
  clumpfile=open(clumpfname,'w')
  #Print some param data and file header
  clumpfile.write('#----------------------------------------------------------------------\n')
  clumpfile.write('# Peak-detection algorithm: Starlink Fellwalker (Berry+2014)\n')
  clumpfile.write('# Params: RMS=%.1f (=sqrt(mean(cts))\n' % (rms))
  if args.frms is not None: 
    clumpfile.write('#         MinHeight=%d (=%.1f*RMS)\n' % (minheight,args.frms))
  else: 
    clumpfile.write('#         MinHeight=%d (=%.2f*max_cts)\n' % (minheight,args.ffrac))
  clumpfile.write('#----------------------------------------------------------------------\n')
  hfmt='#%3s '+6*'%8s '+'%10s '+'\n'
  clumpfile.write(hfmt % ('ID','phi_p','theta_p','phi_c','theta_c','dphi','dtheta','peak_cts'))

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

  #Print out
  fmt='%4d '+6*'%8.3f '+'%10.0f '
  scipy.savetxt(clumpfile,np.array([pid,phipeak,thetapeak,phipeakc,thetapeakc,dphi,dtheta,cheight]).T,fmt=fmt)

  #Plot detected clump peaks
  m.scatter(xpeak,ypeak,c='none',edgecolor='k',s=20,zorder=99)

  #Save current figure 
  fig.savefig(figname)

  #If flag is set, plot new figure indicating pixels associated to each clump
  if args.saveclumps:
    #Convert clump pixel mask to fits
    print 'Converting pixel map from NDF to fits'
    os.system('%s/convert/ndf2fits _zp_cmask _zp_cmask.fits' % (starlink_path)) 

    #Read in fits image
    cmask_dat=pyfits.getdata('_zp_cmask.fits',0)
    inds=np.arange(npix)+1
    xpix_2d,ypix_2d=np.meshgrid(inds,inds)
    xinds,yinds=xpix_2d.flatten(),ypix_2d.flatten()
    pcts_1d=zi.flatten()
    cmask_1d=cmask_dat.flatten()

    #Save only pixels inside the FWHM of the peak and with counts>minheight
    mmask=np.zeros(cmask_1d.size,dtype=bool)
    for kk in range(pid.size):
      #Combine mask
      mmask=mmask | ( (cmask_1d==pid[kk]) & (pcts_1d>=0.5*cheight[kk]))

    #Keep only data inside FWHM
    xcmask,ycmask,phicmask,thetacmask = pix_convert.get_phys_from_pix(xinds[mmask],yinds[mmask])
    cmask_1d=cmask_1d[mmask]
 
    #Plot identified clumps on top of pole count map
    color_cycle=['red','mediumblue','orange','forestgreen','darkorchid','cyan','dodgerblue','teal','magenta','deeppink','maroon']
    for kk in range(pid.size):
      pmask=(cmask_1d==pid[kk])
      m.plot(xcmask[pmask],ycmask[pmask],color=color_cycle[kk],mec='None',ms=5,marker='o',alpha=0.3)

    #Print  pixel data to output pole list file
    head='%s %10s %10s' % ('IDpole','phi_pole','theta_pole')
    sorted_mask=cmask_1d.argsort()
    scipy.savetxt(clumppixfname,np.array([cmask_1d[sorted_mask],phicmask[sorted_mask],thetacmask[sorted_mask]]).T,fmt='%8d %10.4f %10.4f',header=head)

    cfigname='%s.%s.%s.%s.pls.%s' % (figname_root,mode,proj[:3],pmode,args.fig)
    fig.savefig(cfigname)

    #Remove auxiliary files
    os.system('rm -f _zp*')

  if args.show: plt.show()
  else: fig.clf()
