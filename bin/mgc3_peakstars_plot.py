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
import myutils
import mgc3_lib
import gzip
#----------------------------------------------
__version__ = '1.0'
__docformat__ = "reredtext en"
#
parser = argparse.ArgumentParser(description='Plot stars associated to each of the poles in a list')
parser.add_argument('parfile',metavar='parameter_file',help='Input catalogue parameter file',action='store',nargs=1)
parser.add_argument('infile',metavar='inpstfile',help='Input file containing catalogue for stars associated to peaks',nargs=1,action='store')
parser.add_argument("-l", "--llist", action="store_true",help='Take inpstfile as list of mgc3.cts files')
parser.add_argument('-cat','--catfile',metavar='catfile',help='Full input catalogue. Must be a list when using -l',action='store',default=False,nargs=1)
parser.add_argument('-ext',metavar='outext',help='Output suffix [optional]',action='store',default=['',],nargs=1)
parser.add_argument('-lon0',help='Longitude for Y-axis. Default is 0.', action='store',default=0.,type=np.float)
parser.add_argument('-dlat',help='Spacing between parallels. Default is 20.', action='store',default=30.,type=np.float)
parser.add_argument('-dlon',help='Spacing between meridians. Default is 30.', action='store',default=30.,type=np.float)
parser.add_argument('-xlim',metavar='xo xf',help='Set X limits (space-separated)',action='store',nargs=2,type=np.float)
parser.add_argument('-ylim',metavar='yo yf',help='Set Y limits (space-separated)',action='store',nargs=2,type=np.float)
parser.add_argument('-zlim',metavar='zo zf',help='Set Z limits (space-separated)',action='store',nargs=2,type=np.float)
parser.add_argument('-title',help='Plot title', action='store',default=None)
parser.add_argument('-helio',help='Use heliocentric coords in Aitoff plot', action='store_true',default=False)
parser.add_argument('-f','--fig',help='Output plot type png/eps. Default is png', action='store',default='png',choices=['png','eps','pdf'])
parser.add_argument('-ms',help='Marker size for peak stars. Use ms=0 for fullcat only.',action='store',default=3,type=np.float)
parser.add_argument('-s','--show',help='Show plot in window. Default is False', action='store_true',default=False)
parser.add_argument('-ic','--idcol',help='Column containing stream ID (counting from 1). Default is last col.', action='store',default=-1,type=np.int)


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

#If optional input catfile given, check consistency with infile
if args.catfile:
  if not args.llist:
     cat_list=[args.catfile[0],]
  else:
     cat_list=scipy.genfromtxt(args.catfile[0],dtype='S')
     if len(cat_list)!=len(file_list):
        sys.exit('WARNING: CAT_LIST(len=%d) must have the same length as INFILE_LIST(len=%d).\nExiting...' 
                  % (len(cat_list),len(file_list)))
     if np.ndim(cat_list)==0: cat_list=array([cat_list,])

parfile = args.parfile[0]
print 'Reading Parameter file %s ...' % (parfile)
spars=mgc3_lib.parse_pars(parfile)

#Parse xyz-limits if given
if args.xlim: print 'X limits:', args.xlim
if args.ylim: print 'Y limits:', args.ylim
if args.zlim: print 'Z limits:', args.zlim

for ff in range(len(file_list)):

  infilen=file_list[ff]
  if '.gz' in infilen: infile=gzip.open(infilen)
  else: infile=open(infilen)
  dat=scipy.genfromtxt(infile,comments='#')

  print 'Plotting stars from file: ',infilen

  if spars['par_muas']: fp=1.
  else: fp=1000.
  if spars['pm_muas']: fm=1.
  else: fm=1000.
  unit='(kpc)'

  if args.catfile:
    if 'gz' in cat_list[ff]:
      cat_listf=gzip.open(cat_list[ff],'r')
    else:
      cat_listf=open(cat_list[ff],'r')
    cdat=scipy.genfromtxt(cat_listf,comments='#')
    #Apply whatever filters may be defined by the AUX pars
    mask=np.ones(cdat[:,0].size,dtype='bool')
    for NAUX in range(1,spars['NAUX']+1,1):
      mykey_col='AUX%d_col' % (NAUX)
      mykey_valo='AUX%d_o' % (NAUX)
      mykey_valf='AUX%d_f' % (NAUX)
      #Skip if col=998   
      if spars[mykey_col]!=998:
       print ' Cutting input catalogue with %.1f<%s[%d]<%.1f' % (spars[mykey_valo],
              mykey_col,spars[mykey_col]+1,spars[mykey_valf])
       #Create mask 
       mask_i = (cdat[:,spars[mykey_col]]>spars[mykey_valo]) & (cdat[:,spars[mykey_col]]<spars[mykey_valf])
       #Combine masks
       mask = mask & mask_i
    #Filter
    cdat=cdat[mask,:]
    #Parse
    cl,cb,cparallax=cdat[:,spars['lon_col']],cdat[:,spars['lat_col']],cdat[:,spars['par_col']]
    cmulstar,cmub,cvrad=cdat[:,spars['pm_lon_col']],cdat[:,spars['pm_lat_col']],cdat[:,spars['vrad_col']]
    css=myutils.helio_obj(cl,cb,fp*cparallax,fm*cmulstar,fm*cmub,cvrad,degree=spars['deg'],flag_mulstar=spars['pm_lon_red'])

  figname_root=infilen.replace('.pst','')
  fig1name='%s.xyz%s.%s' % (figname_root,args.ext[0],args.fig)
  fig2name='%s.sph%s.%s' % (figname_root,args.ext[0],args.fig)

  IDpole=dat[:,args.idcol-1]
  l,b,parallax=dat[:,spars['lon_col']],dat[:,spars['lat_col']],dat[:,spars['par_col']]
  mulstar,mub,vrad=dat[:,spars['pm_lon_col']],dat[:,spars['pm_lat_col']],dat[:,spars['vrad_col']]
  
  ss=myutils.helio_obj(l,b,fp*parallax,fm*mulstar,fm*mub,vrad,degree=spars['deg'],flag_mulstar=spars['pm_lon_red'])

  #ss.x,ss.y,ss.z=dat[:,8-1],dat[:,9-1],dat[:,10-1] #for tests only

  npoles=np.int(np.max(IDpole))
 # npoles_id=IDpole.unique().size
  
  print 'Npeaks=',npoles
  cmapp=plt.cm.gist_ncar_r(np.linspace(0.1, 0.9, npoles))  #Upper limit is 0.85 to avoid last colors of the colormap
#  cmapp=plt.cm.spectral(np.linspace(0.1, 0.9, npoles))  #Upper limit is 0.85 to avoid last colors of the colormap
  if npoles<=10:
     cmapp=['darkviolet','orange','lime','royalblue','orchid','red','gray','pink','limegreen','navy']

#     cmapp=['darkviolet','slateblue','deeppink','royalblue','orchid','red','gray','pink','limegreen','navy']
  #cmapp=['darkviolet','orange','lime','royalblue','orchid','red','gray','pink','limegreen','navy']
#     cmapp=['orchid','red','mediumblue','orange','red','royalblue','gray','pink','limegreen','navy']


  fig1=plt.figure(1,figsize=(13,6))
  fig1.subplots_adjust(wspace=0.2,left=0.08,right=0.97)
  ax1=fig1.add_subplot(121)
  ax2=fig1.add_subplot(122)

  fig2=plt.figure(2,figsize=(10,6))
  fig2.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)
  ax4=fig2.add_subplot(111)
  m = Basemap(projection='moll',ax=ax4,lon_0=args.lon0)
  mers=np.arange(0.,360.,args.dlon)
  parls=np.arange(-90.,90.,args.dlat)
  m.drawmeridians(mers,color='gray')
  m.drawparallels(parls,color='gray')
  m.drawmapboundary()

  c_props={'ms':0.5,'zorder':0,'alpha':0.7}
  s_props={'ms':args.ms,'zorder':1,'alpha':1.,'mec':'None'}

  for kk in range(npoles):
   mask= (IDpole==kk+1)
   #---Cartesian---------
   ax1.plot(-ss.x[mask],ss.y[mask],'.',color=cmapp[kk],**s_props)
   #---
   ax2.plot(-ss.x[mask],ss.z[mask],'.',color=cmapp[kk],**s_props)
   #---Aitoff--------
   if args.helio: xmoll,ymoll=m(ss.l[mask],ss.b[mask])
   else:          xmoll,ymoll=m(ss.phi[mask],ss.theta[mask])
   m.plot(xmoll,ymoll,'.',color=cmapp[kk],**s_props)

  if args.catfile:
   ax1.plot(-css.x,css.y,'k.',**c_props)
   ax2.plot(-css.x,css.z,'k.',**c_props)
   cxmoll,cymoll=m(css.phi,css.theta)
   m.plot(cxmoll,cymoll,'k.',**c_props)

  #Aitoff plot labels
  lmers=mers[mers!=args.lon0+180.]
  bmers=0.+np.zeros_like(lmers)
  bpars=parls[(np.abs(parls)!=90.) & (np.abs(parls)!=0.)]
  lpars=args.lon0+180.+np.zeros_like(bpars)
  xl,yl=m(lmers,bmers)
  for ii in range(lmers.size): ax4.text(xl[ii],yl[ii],"%.0f$^o$" % (lmers[ii]),horizontalalignment='center',
                                         fontsize=11,color='gray')
  xl,yl=m(lpars,bpars)
  for ii in range(bpars.size): ax4.text(xl[ii],yl[ii],"  %+.0f$^o$" % (bpars[ii]),verticalalignment='center',
                                           horizontalalignment='left',fontsize=11,color='gray')
  #If flags are set, force xyz-limits
  if args.xlim: 
   ax1.set_xlim(args.xlim)
   ax2.set_xlim(args.xlim)
  if args.ylim: 
   ax1.set_ylim(args.ylim)
  if args.zlim: 
   ax2.set_ylim(args.zlim)

  ax1.set_xlabel('X '+unit)
  ax2.set_xlabel('X '+unit)
  ax1.set_ylabel('Y '+unit)
  ax2.set_ylabel('Z '+unit)
  if args.title: ax1.set_title(args.title)
  if args.title: ax4.set_title(args.title)
  fig1.savefig(fig1name)
  fig2.savefig(fig2name)

  if args.show: plt.show()
  else: 
    fig1.clf()
    fig2.clf()
