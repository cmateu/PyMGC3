#!/usr/bin/env python
import sys
import scipy
import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#############################################################################
#Copyright (c) 2013 - 2014, Cecilia Mateu
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without 
#modification, are permitted provided that the following conditions are met:
#
#   Redistributions of source code must retain the above copyright notice, 
#      this list of conditions and the following disclaimer.
#   Redistributions in binary form must reproduce the above copyright notice, 
#      this list of conditions and the following disclaimer in the 
#      documentation and/or other materials provided with the distribution.
#   The name of the author may not be used to endorse or promote products 
#      derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
#BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
#OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
#AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
#WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#POSSIBILITY OF SUCH DAMAGE.
#############################################################################

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

