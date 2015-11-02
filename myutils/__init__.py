#!/usr/bin/env python
import sys
import scipy
import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import bovy_coords as bovyc
import gzip 
import matplotlib.colors 
import matplotlib.cm 

#############################################################################
# VERSION: 
# 2015/02 - v2 
#           Obj Class added for plotting peak stars utility
#############################################################################

def get_header_line(filename):

 #Open file
 if '.gz' in filename:
     ifile=gzip.open(filename,'r')
     filename=filename.replace('.gz','')
 else: ifile=open(filename,'r')

 line='0'
 head=[]
 while len(line)>0:
  #Read lines one by one and chomp (i.e. remove new line character at the end of line)
  line=ifile.readline().rstrip('\r\n')
  if '#' in line: head.append(line)
  else: break

 if len(head)==0:
   head=['# ',]

 return head

def aitoff_plot(l,b,ax=None,show=True,mt='k.',ms=1.):

  if ax is None:
    fig = plt.figure(figsize=(7.4,7.5))
    ax = fig.add_subplot(111, projection='aitoff')

  ll=l
  ll[ll>=180.]=ll[ll>180.]-360.
  ll=np.radians(ll)
  bb=np.radians(b)

  ax.plot(ll,bb,mt,ms=ms)
  if show: plt.show()
  return ax

def spherical_plot(lon,lat,deg=True,c='k',ms=3.,ax=None,show=True):

  if ax is None:
   fig = plt.figure(figsize=(7.4,7.5))
   ax = fig.add_subplot(111, projection='3d')
  
  u = np.linspace(0, 2*np.pi, 100)
  v = np.linspace(0, np.pi, 100)
  
  x = 10 * np.outer(np.cos(u), np.sin(v))
  y = 10 * np.outer(np.sin(u), np.sin(v))
  z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))
  ax.plot_wireframe(x, y, z, rstride=10, cstride=10, color='gray', lw=1.)
  
  if deg: lon,lat=np.radians(lon),np.radians(lat)
  xx = 10*np.cos(lon)*np.cos(lat)
  yy = 10*np.sin(lon)*np.cos(lat)
  zz = 10*np.sin(lat)
  ax.scatter(xx,yy,zz,c=c,s=ms)

  if show: plt.show()

  return ax

def uniform_spherical_data(N):
 xx=np.random.random(N)
 lon=2*np.pi*np.random.random(N)
 lat=np.arcsin(2*xx-1)
 return (lon,lat)

#========================================================================

class helio_obj:
  #Default values for Sun's position and velocity are from Brown et al. 2005
  def __init__(self,l,b,parallax,mulstar,mub,vrad,flag_mulstar=True,xyz_sun=[-8.5,0.,0.],vel_sun=[10.3,232.6,5.9],verbose=False,degree=False):

    #Degree2Radian conversion
    if degree: _d2r=np.pi/180.
    else: _d2r=1.

    #Save inputs
    self.l,self.b,self.parallax=l,b,parallax
    self.mub,self.vrad,self.Rhel=mub,vrad,1000./self.parallax  #if par in mas/muas, rhel in pc/kpc
    if flag_mulstar:
       self.mulstar=mulstar
       self.mul=mulstar/np.cos(self.b*_d2r)
    else:
       self.mul=mulstar
       self.mulstar=self.mul*np.cos(self.b*_d2r)

    #Bovy's library assumes Sun's position is positive. Flip X-axis if xsun<0
    if xyz_sun[0]<0.: sign=-1
    else: sign=+1

    #Save Sun's position
    self.xsun, self.ysun, self.zsun= xyz_sun
    self.vxsun, self.vysun, self.vzsun= vel_sun
    self.xsun, self.vxsun = sign*self.xsun, sign*self.vxsun

    #Convert to heliocentric cartesian. Bovy's library assumes Sun's position is positive
    #tuple output, no .T needed
    if verbose: print 'Converting Heliocentric Galactic Spherical to Heliocentric Cartesian coords...'
    m=bovyc.lbd_to_XYZ(self.l,self.b,self.Rhel,degree=degree)
    self.xhel,self.yhel,self.zhel=m.T
    m=bovyc.vrpmllpmbb_to_vxvyvz(self.vrad,self.mulstar,self.mub,self.l,self.b,self.Rhel,XYZ=False,degree=degree)
    self.vxhel,self.vyhel,self.vzhel=m.T
    #m=bovyc.galcenrect_to_XYZ(self.x,self.y,self.z,Xsun=self.xsun,Ysun=self.ysun,Zsun=self.zsun)a

    #Convert Heliocentric Cartesian to Galactocentric Cartesian
    if verbose: print 'Converting Heliocentric Cartesian to Galactocentric Cartesian coords...'
    m=bovyc.XYZ_to_galcenrect(self.xhel,self.yhel,self.zhel,Xsun=self.xsun,Ysun=self.ysun,Zsun=self.zsun)
    self.x,self.y,self.z=m
    m=bovyc.vxvyvz_to_galcenrect(self.vxhel,self.vyhel,self.vzhel,vsun=[self.vxsun, self.vysun, self.vzsun])
    self.vx,self.vy,self.vz=m

    #Compute Galactocentric Spherical
    if verbose: print 'Converting Galactocentric Cartesian to Spherical coords...'
    m=bovyc.XYZ_to_lbd(self.x,self.y,self.z,degree=degree)
    self.phi,self.theta,self.Rgal=m.T
    self.phi=(self.phi+180.) % 360.    

#--------------color maps---------------------------------------------------

def create_cmap_from_list(colorlist,cmap_name):

  cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list(cmap_name,colorlist)
  matplotlib.cm.register_cmap(cmap=cmap2)
  return cmap2


#Implementations from Paul Tol's SRON notes on Colour Maps

def get_sron_rainbow(N=12):

 if N==12:
   colorlist=['#781C81','#413B93','#4065B1','#488BC2','#55A1B1','#63AD99',
   #colorlist=['#413B93','#4065B1','#488BC2','#55A1B1','#63AD99',
              '#7FB972','#B5BD4C','#D9AD3C','#E68E34','#E6642C','#D92120']
   cmap=create_cmap_from_list(colorlist,'sron_rainbow12')

 return cmap



