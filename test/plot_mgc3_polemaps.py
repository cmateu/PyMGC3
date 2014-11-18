#!/usr/bin/env python
import pylab as plt
import scipy
import numpy as np
from mpl_toolkits.basemap import Basemap
import sys
import argparse
import scipy.ndimage

parser = argparse.ArgumentParser()
parser.add_argument('infile',metavar='infile',help='Input ascii file to be converted to fits and compressed with fpack',nargs=1,action='store')
parser.add_argument("-l", "--llist", action="store_true",help='Take infile as list of mgc3.cts files')
parser.add_argument('-m',help='Plot mGC3/nGC3/GC3 pole count map. Default is mGC3', action='store',default='mGC3',choices=['mGC3','nGC3','GC3'])
parser.add_argument('-f','--fig',help='Output plot type png/eps. Default is png', action='store',default='png',choices=['png','eps'])
parser.add_argument('-proj',help='Projection npaeqd/ortho/mollweide. Default is npaeqd', action='store',default='npaeqd',choices=['npaeqd','ortho','mollweide'])
parser.add_argument('-c','--contour',help='Plot pole-count contour map instead of raw grid.', action='store_true',default=False)
parser.add_argument('-t','--twohemispheres',help='Plot both hemispheres in pole-count map.', action='store_true',default=False)

#---------Parse----------------------------
args = parser.parse_args()

#Parse inputs
if not args.llist:
 print 'Reading file: %s' % (args.infile)
 file_list=[args.infile[0],]
else:
 print 'Reading input files from list file: ', file_list
 file_list=genfromtxt(args.infile[0],dtype='S')
 if ndim(file_list)==0: file_list=array([file_list,])

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

ori='vertical'
ori='horizontal'
ni=0

colormap=plt.cm.jet
for infilen in file_list:

  phio,thetao,pole_ctso=pdat=scipy.genfromtxt(infilen,comments='#',usecols=(0,1,counts_col),unpack=True)
  figname_root=infilen.replace('.mgc3.cts','')
  figname='%s.%s.%s.%s' % (figname_root,mode,pmode,args.fig)
  print 'Output filename:', figname
  sys.exit()

  phi=np.append(phio,(phio+180.) % 360.)
  theta=np.append(thetao,-thetao)
  pole_cts=np.append(pole_ctso,pole_ctso) 

  #----------------Pole count map------------------------------
  nx,ny=2,3

  mer_grid=[0.,360.,20.]
  par_grid=[-90.,+90.,30.]
  ms=50.
  proj='npaeqd'
 
  flag_raw=True
  fig=plt.figure(1,figsize=(14,8))
  if flag_raw: 
     ax=fig.add_subplot(nx,ny,1)
     m = Basemap(projection=proj,boundinglat=0.,lon_0=0,resolution='l',ax=ax)
     m.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]))
     m.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]))
     x,y=m(phi,theta)
     c=m.scatter(x,y,c=pole_cts,edgecolor='none',s=ms,cmap=colormap)
     ax.set_title('%s (RAW)' % (mode_ori))
     #plt.colorbar(c,pad=0,orientation=ori,format='%d')

  npix=250
  xi = np.linspace(np.min(x),np.max(x),npix)
  yi = np.linspace(np.min(y),np.max(y),npix)
  zi = plt.griddata(x,y,pole_cts,xi,yi) #,'nn')

  #my_extent=[xi[0],xi[-1],yi[0],yi[-1]]
  clevels=30
  ax=fig.add_subplot(nx,ny,2)
  m = Basemap(projection=proj,boundinglat=-10,lon_0=0,resolution='l',ax=ax)
  m.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]))
  m.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]))
  m.contourf(xi,yi,zi,clevels,cmap=colormap)
  ax.set_title('%s (raw contours)' % (mode_ori))
  #ax.imshow(zi,interpolation='nearest',origin='lower',aspect='auto',cmap=colormap) #no funciona con basemap todavia

  #--------------------UNSHARP MASKING--------------------------------------------------------------------
  ax=fig.add_subplot(nx,ny,3)
  m = Basemap(projection='npaeqd',boundinglat=-10,lon_0=0,resolution='l',ax=ax)
  m.drawmeridians(np.arange(0.,360.,20.))
  m.drawparallels(np.arange(-90.,90.,30.))
  nsm=30
  zi_smooth=scipy.ndimage.median_filter(zi,size=(nsm,nsm),mode='wrap')
  m.contourf(xi,yi,zi_smooth,25,cmap=colormap)

  ax=fig.add_subplot(nx,ny,4)
  zi_sharp=zi-zi_smooth
  phii,thetai=m(xi,yi)
  ax.imshow(zi_sharp,interpolation='nearest',origin='lower',aspect='auto',cmap=colormap)

  #Reshape matrix into 1D array
  M=zi_sharp
  M_1d=M.flatten()
  xs=np.reshape(np.array(len(yi)*list(xi)),np.shape(M))
  ys=np.reshape(np.array(len(xi)*list(yi)),np.shape(M.T))
  ys=ys.transpose()
  xs=xs.flatten()
  ys=ys.flatten()
  phis,thetas=m(xs,ys,inverse=True)

  ax=fig.add_subplot(nx,ny,5)
  m = Basemap(projection='npaeqd',boundinglat=0,lon_0=0,resolution='l',ax=ax)
  m.drawmeridians(np.arange(0.,360.,20.))
  m.drawparallels(np.arange(-90.,90.,30.))
  my_extent=[xi[0],xi[-1],yi[0],yi[-1]]
  #ax.imshow(zi_sharp,interpolation='nearest',origin='lower',aspect='auto',cmap=colormap,extent=my_extent)
  #m.contourf(xi,yi,zi_sharp,30,cmap=colormap)
  c=m.scatter(xs,ys,c=M_1d,edgecolor='none',s=50,cmap=colormap)


  proj='ortho'
  print proj
  ax=fig.add_subplot(nx,ny,6)
  m = Basemap(projection=proj,lat_0=50,lon_0=0.,resolution='l',ax=ax,area_thresh = 1000.)
  m.drawmeridians(np.arange(0, 360, 20))
  m.drawparallels(np.arange(-90, 90, 20))
  m.drawmapboundary()
  xso,yso=m(phis,thetas)
  c=m.scatter(xso,yso,c=M_1d,edgecolor='none',s=20,cmap=colormap)


  fig.savefig(figname_a)
  plt.show()
  fig.clf()
