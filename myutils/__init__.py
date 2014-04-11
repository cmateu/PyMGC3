#!/usr/bin/env python
import sys
import scipy
import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def get_header_line(filename):

 file=open(filename,'r')
 line='0'
 head=[]
 while len(line)>0:
  #Read lines one by one and chomp (i.e. remove new line character at the end of line)
  line=file.readline().rstrip('\r\n')
  if '#' in line: head.append(line)
  else: break

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

