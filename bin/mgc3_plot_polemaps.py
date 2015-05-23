#!/usr/bin/env python
import pylab as plt
import scipy
import numpy as np
from mpl_toolkits.basemap import Basemap
import sys
import argparse
import scipy.ndimage

__version__ = '1.1'
__docformat__ = "reredtext en"
__what__= sys.argv[0]+": This program detects peaks in pole count maps using the Fellwalker algorithm (starlink implementation)"
#
parser = argparse.ArgumentParser(description='Plot mGC3/nGC3/GC3 pole count maps')
parser.add_argument('infile',metavar='infile',help='Input file containing pole count maps (*.cts file)',nargs=1,action='store')
parser.add_argument("-l", "--llist", action="store_true",help='Take infile as list of mgc3.cts files')
parser.add_argument('-m',help='Plot mGC3/nGC3/GC3/mGC3hel pole count map. Default is mGC3', action='store',default='mGC3',choices=['mGC3','nGC3','GC3','mGC3hel'])
parser.add_argument('-f','--fig',help='Output plot type png/eps. Default is png', action='store',default='png',choices=['png','eps','pdf'])
parser.add_argument('-proj',help='Projection npaeqd/ortho/mollweide. Default is npaeqd', action='store',default='npaeqd',choices=['npaeqd','ortho','moll'])
parser.add_argument('-log',help='Plot pole-count map in log-scale', action='store_true',default=False)
parser.add_argument('-lon0',help='Longitude for Y-axis. Default is 0.', action='store',default=0.,type=np.float)
parser.add_argument('-lat0',help='Bounding latitude for plot. Default is 90.', action='store',default=0.,type=np.float)
parser.add_argument('-dlat',help='Spacing between parallels. Default is 20.', action='store',default=20.,type=np.float)
parser.add_argument('-dlon',help='Spacing between meridians. Default is 20.', action='store',default=30.,type=np.float)
parser.add_argument('-latmax',help='Max latitude upto which meridians are drawn. Default is 80.', action='store',default=80.,type=np.float)
parser.add_argument('-vmin',help='Min counts for color-scale. Default is min(cts)', action='store',default=None,type=np.float)
parser.add_argument('-vmax',help='Max counts for color-scale. Default is max(cts)', action='store',default=None,type=np.float)
parser.add_argument('-ms',help='Marker size. Default: 15/40 for npaeqd/ortho.', action='store',default=-1.,type=np.float)
parser.add_argument('-c','--contour',help='Plot pole-count contour map instead of raw grid.', action='store_true',default=False)
parser.add_argument('-t','--twohemispheres',help='Plot both hemispheres in pole-count map.', action='store_true',default=False)
parser.add_argument('-s','--show',help='Show plot in window. Default is False', action='store_true',default=False)
parser.add_argument('-title',help='Plot title', action='store',default=None)
parser.add_argument('-pls',metavar='PLSFILE',help='Overplot poles from peakdetect output file (.pls)', action='store',default=None)
parser.add_argument('-al','--alpha',help='Clump transparency. Default 0.4', action='store',default=0.4,type=np.float)
parser.add_argument('-bw',help='Use grayscale colormap to plot PCMs. Default False (uses jet colormap)', action='store_true',default=False)
parser.add_argument('-ext',metavar='outfile_ext',help='Output suffix [optional]. If given output will be infile.outfile_ext.mgc3.pst',action='store',default=['',],nargs=1)

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
if 'hel' in mode:    counts_col=4-1
elif 'mgc3' in mode: counts_col=3-1
elif 'ngc3' in mode: counts_col=6-1
elif 'gc3'  in mode: counts_col=5-1

#Parse raw/contour mode-------------------------
if args.contour: 
  pmode='c'
  print 'Plotting contour pole-count map'
else: 
  pmode='r'
  print 'Plotting raw pole-count map'
if args.log:
  print 'Plotting pole-count map in log-scale'
  pmode=pmode+'l'

if args.twohemispheres:
  print 'Plotting both hemispheres in pole-count map'
else: print 'Plotting one hemisphere in pole-count map.'

print 'Plotting using projection:', args.proj

ori='vertical'
ori='horizontal'
ni=0

colormap=plt.cm.jet
if args.bw: colormap=plt.cm.gray

for infilen in file_list:

  #Default title----------------------------------
  args.title=infilen

  phio,thetao,pole_ctso=pdat=scipy.genfromtxt(infilen,comments='#',usecols=(0,1,counts_col),unpack=True)
  figname_root=infilen.replace('.mgc3.cts',args.ext[0])  #works well if args.ext is empty
  figname='%s.%s.%s.%s.%s' % (figname_root,mode,args.proj[:3],pmode,args.fig)
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
  nx,ny=1,1

  mer_grid=[0.,360.,args.dlon]
  #par_grid=[-90.,+90.,args.dlat]
  par_grid=[-args.dlat,+90.,args.dlat]

  if 'npa' in args.proj or 'moll' in args.proj:
    #For npa and moll projections, plot map as viewed from lon0 only
    fig=plt.figure(1,figsize=(8,8))
    dw=0.8
    wo=(1.-dw)/2.
    wyo=0.05*wo
    fig.subplots_adjust(left=wo,right=dw+wo,top=dw+wo,bottom=wyo)
    nrow,ncol,nplot=1,1,1
    opts=[(1,args.lon0),] 
    proj_dict={'boundinglat':args.lat0,'resolution':'l'}
    if args.ms==-1: ms=70.
    else: ms=args.ms
  else:
    #For ortho projection, plot map as viewed from lon=0 and lon0+180
    fig=plt.figure(1,figsize=(12,6))
    nrow,ncol=1,2
    opts=[(1,args.lon0),(2,args.lon0+180.)] 
    proj_dict={'boundinglat':args.lat0,'resolution':'l','lat_0':50.,'area_thresh':1000.}
    if args.ms==-1: ms=50.
    else: ms=args.ms

  for ii,l0 in opts:
    ax=fig.add_subplot(nrow,ncol,ii)
    m = Basemap(projection=args.proj,lon_0=l0,ax=ax,**proj_dict)
    m.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]),color='lightgrey',lw=2.)
    m.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]),color='lightgrey',lw=2.)
    m.drawmapboundary()

    x,y=m(phi,theta)

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
       npix=300
       clevels=30
       xi = np.linspace(np.min(x),np.max(x),npix)
       yi = np.linspace(np.min(y),np.max(y),npix)
       zi = plt.griddata(x,y,pole_cts,xi,yi) #,'nn')
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
    #Labels and such
    if 'npa' not in args.proj: ax.set_title('%s pole-counts' % (mode_ori))

  #If given, read in pls file
  if args.pls is not None:
   poleIDs,phis,thetas=scipy.genfromtxt(args.pls,unpack=True,usecols=(0,1,2))
   u_pid=np.unique(poleIDs)
   #Define colormap consistently with peakdetect
   cmapp=plt.cm.gist_ncar_r(np.linspace(0.1, 0.9, u_pid.size))  #Upper limit is 0.85 to avoid last colors of the colormap
   if u_pid.size<=10:
     cmapp=['darkviolet','orange','lime','royalblue','orchid','red','gray','pink','limegreen','navy']
   for kk in np.arange(u_pid.size):
     #Project pole coords
     xpoles,ypoles=m(phis,thetas)
     #plot current peak only
     idmask=poleIDs==u_pid[kk]
     m.scatter(xpoles[idmask],ypoles[idmask],c=cmapp[kk],edgecolors='none',s=20,marker='o',alpha=args.alpha)
     #Peak ID labels
     u_xpeak,u_ypeak=np.median(xpoles[idmask]),np.median(ypoles[idmask])
     m.scatter(u_xpeak,u_ypeak,c='w',alpha=0.5,edgecolor='k',s=110,zorder=100)
     ax.text(u_xpeak,u_ypeak,'%d' % (u_pid[kk]),fontsize=7,color='black',
                horizontalalignment='center',verticalalignment='center',zorder=101)


  #Plot colorbar
  if 'npa' in args.proj:
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

  fig.savefig(figname)
  if args.show: plt.show()
  else: fig.clf()


