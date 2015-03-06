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
parser.add_argument('-c','--contour',help='Plot pole-count contour map instead of raw grid.', action='store_true',default=False)
parser.add_argument('-t','--twohemispheres',help='Plot both hemispheres in pole-count map.', action='store_true',default=False)
parser.add_argument('-s','--show',help='Show plot in window. Default is False', action='store_true',default=False)
parser.add_argument('-sc','--saveclumps',help='Plot and save poles associated to each peak.', action='store_true',default=False)
peakargs = parser.add_mutually_exclusive_group()
peakargs.add_argument('-fr','--frms',help='If set, min peak height is frms*RMS', action='store',type=np.float)
peakargs.add_argument('-ff','--ffrac',help='Default option. Min peak height is fmax*max_pole_counts. Default fmax=0.6', action='store',default=0.6,type=np.float)
parser.add_argument('-mj','--maxjump',help='Fellwalker MaxJump param, neighbourhood radius to search for +gradient. Default 20.', action='store',default=20,type=np.float)
parser.add_argument('-al','--alpha',help='Clump transparency. Default 0.3', action='store',default=0.3,type=np.float)
parser.add_argument('-fx','--fwxm',help='Store pixels with cts>fwxm*maxpeak. Default is 0.5 (=FWHM)', action='store',default=0.5,type=np.float)


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
  #If log-flag is set, do everything with log(counts)
  if args.log: 
     print 'Do peak detection on linear scale but show plot on log scale'
     #pole_ctso=np.log10(pole_ctso)
     pmode=pmode+'l'

  #Output figure and file names
  figname_root=infilen.replace('.mgc3.cts',args.ext[0])  #works well if args.ext is empty
  figname='%s.%s.%s.%s.%s' % (figname_root,mode,proj[:3],pmode,args.fig)
  clumpfname='%s.%s.peak.pls' % (figname_root,mode)
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
  par_grid=[-args.dlat,+90.,args.dlat]

  #For npa and moll projections, plot map as viewed from lon0 only
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
  #
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
  os.system('%s/cupid/findclumps in=_zpndf.sdf out=_zp_cmask method=fellwalker outcat=_zp_clumps rms=%f config="fellwalker.maxjump=%.0f" ' % (starlink_path,rms,args.maxjump)) 

  #Read-in output table with identified clumps
  clumpdat=pyfits.open('_zp_clumps.FIT')[1].data                    
  xpix,ypix,cheight=clumpdat.field('Peak1'),clumpdat.field('Peak2'),clumpdat.field('Peak')
  xpix_c,ypix_c,dx_pix,dy_pix=clumpdat.field('Cen1'),clumpdat.field('Cen2'),clumpdat.field('Size1'),clumpdat.field('Size2')
  pid=np.arange(xpix.size) + 1  #Peak IDs, must start from 1

  #Keep only peaks with counts > minheight and with sizes>=1 (clumps have to be at least 1 pixel in size)
  mask=(cheight>=minheight) & (dx_pix>=1.) & (dy_pix>=1.)
  if not mask.any(): 
    print 'No peaks above minheight or with size>=1pix found'
    continue
  xpix,ypix,cheight=xpix[mask],ypix[mask],cheight[mask]
  xpix_c,ypix_c,dx_pix,dy_pix,pid=xpix_c[mask],ypix_c[mask],dx_pix[mask],dy_pix[mask],pid[mask]
  newid=np.arange(xpix.size) + 1  #Rename so IDs will be consecutive numbers starting from 1

  #print pars on screen
  print '#----------------------------------------------------------------'
  print '# Peak-detection algorithm: Starlink Fellwalker (Berry 2014)'
  print '# Params: RMS=%.1f (=sqrt(mean(cts))' % (rms)
  print '#         FellWalker.MaxJump=%.0f'
  if args.frms is not None:
    print '#         MinHeight=%d (=%.1f*RMS)' % (minheight,args.frms)
  else:
    print '#         MinHeight=%d (=%.2f*max_cts)' % (minheight,args.ffrac)
  print '#         FWXM=%.2f (stored pixels associated to each clump)' % (args.fwxm)
  print '#----------------------------------------------------------------'
  print '# NC=%d clumps saved ' % (xpix.size)
  for kk in range(pid.size): print '# Clump oldID=%3d, newID=%3d, Ncts=%d' % (pid[kk],newid[kk],cheight[kk])
  print '#----------------------------------------------------------------'

  #Open output file to store clump coords
  clumpfile=open(clumpfname,'w')
  #Print some param data and file header
  clumpfile.write('#----------------------------------------------------------------------\n')
  clumpfile.write('# Peak-detection algorithm: Starlink Fellwalker (Berry 2014)\n')
  clumpfile.write('# Params: RMS=%.1f (=sqrt(mean(cts))\n' % (rms))
  clumpfile.write('#         FellWalker.MaxJump=%.0f\n' % (args.maxjump))
  if args.frms is not None: 
    clumpfile.write('#         MinHeight=%d (=%.1f*RMS)\n' % (minheight,args.frms))
  else: 
    clumpfile.write('#         MinHeight=%d (=%.2f*max_cts)\n' % (minheight,args.ffrac))
  clumpfile.write('#         FWXM=%.2f (stored pixels associated to each clump)\n' % (args.fwxm))
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

    xcmask,ycmask,phicmask,thetacmask = pix_convert.get_phys_from_pix(xinds,yinds)
 
    #Plot identified clumps on top of pole count map and print out
    cmapp=plt.cm.gist_ncar(np.linspace(0., 0.9, pid.size))  #Upper limit is 0.85 to avoid last colors of the colormap
    cmapp=plt.cm.gist_ncar_r(np.linspace(0.1, 0.9, pid.size))  #Upper limit is 0.85 to avoid last colors of the colormap
    if pid.size<10:
     cmapp=['mediumblue','orange','lime','orchid','red','royalblue','gray','pink','limegreen','navy']
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
