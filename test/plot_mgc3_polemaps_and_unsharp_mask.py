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

args = parser.parse_args()
mode=args.m.lower()

#Parse inputs
if not args.llist:
 print 'Reading: %s' % (args.infile)
 file_list=[args.infile[0],]
else:
 print 'Reading input files from list file: ', file_list
 file_list=genfromtxt(args.infile[0],dtype='S')
 if ndim(file_list)==0: file_list=array([file_list,])

#Mode
print 'Pole counts plotted: ', mode
if 'mgc3' in mode:   counts_col=3-1
elif 'ngc3' in mode: counts_col=6-1
elif 'gc3'  in mode: counts_col=5-1

ori='vertical'
ori='horizontal'
ni=0

colormap=plt.cm.jet
for infilen in file_list:


  phio,thetao,pole_ctso=pdat=scipy.genfromtxt(infilen,comments='#',usecols=(0,1,counts_col),unpack=True)
  figname=infilen.replace('mgc3.cts',mode)
  figname_a='%s.npaeq.%s' % (figname,args.fig)
  figname_b='%s.ortho.%s' % (figname,args.fig)

  phi=np.append(phio,(phio+180.) % 360.)
  theta=np.append(thetao,-thetao)
  pole_cts=np.append(pole_ctso,pole_ctso) 

  #----------------Pole count map------------------------------
  nx,ny=2,3
  fig=plt.figure(1,figsize=(14,8))
  ax=fig.add_subplot(nx,ny,1)

  m = Basemap(projection='npaeqd',boundinglat=-10,lon_0=0,resolution='l',ax=ax)
  m.drawmeridians(np.arange(0.,360.,20.))
  m.drawparallels(np.arange(-90.,90.,30.))
  x,y=m(phi,theta)
  c=m.scatter(x,y,c=pole_cts,edgecolor='none',s=50,cmap=colormap)
  #plt.colorbar(c,pad=0,orientation=ori,format='%d')

  npix=250
  xi = np.linspace(np.min(x),np.max(x),npix)
  yi = np.linspace(np.min(y),np.max(y),npix)
  zi = plt.griddata(x,y,pole_cts,xi,yi) #,'nn')

  #my_extent=[xi[0],xi[-1],yi[0],yi[-1]]
  ax=fig.add_subplot(nx,ny,2)
  m = Basemap(projection='npaeqd',boundinglat=-10,lon_0=0,resolution='l',ax=ax)
  m.drawmeridians(np.arange(0.,360.,20.))
  m.drawparallels(np.arange(-90.,90.,30.))
  #ax.imshow(zi,interpolation='nearest',origin='lower',aspect='auto',cmap=colormap)
  m.contourf(xi,yi,zi,30,cmap=colormap)

  ax=fig.add_subplot(nx,ny,3)
  m = Basemap(projection='npaeqd',boundinglat=-10,lon_0=0,resolution='l',ax=ax)
  m.drawmeridians(np.arange(0.,360.,20.))
  m.drawparallels(np.arange(-90.,90.,30.))
  nsm=10
  zi_smooth=scipy.ndimage.median_filter(zi,size=(nsm,nsm),mode='wrap')
  m.contourf(xi,yi,zi_smooth,25,cmap=colormap)

  ax=fig.add_subplot(nx,ny,4)
  zi_sharp=zi-zi_smooth
  #ax.imshow(zi_sharp,interpolation='nearest',origin='lower',aspect='auto',cmap=colormap)
  ax.imshow(np.log10(zi_sharp),interpolation='nearest',origin='lower',aspect='auto',cmap=colormap)

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
  m = Basemap(projection='npaeqd',boundinglat=-10,lon_0=0,resolution='l',ax=ax)
  m.drawmeridians(np.arange(0.,360.,20.))
  m.drawparallels(np.arange(-90.,90.,30.))
  m.contourf(xi,yi,zi_sharp,30,cmap=colormap)
  #c=m.scatter(xs,ys,c=M_1d,edgecolor='none',s=50,cmap=colormap)


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
