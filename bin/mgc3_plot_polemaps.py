#!/usr/bin/env python
import pylab as plt
import scipy
import scipy.interpolate
import numpy as np
from mpl_toolkits.basemap import Basemap
import sys
import argparse
import scipy.ndimage
import myutils
import myutils.newpy_colormaps as newcmap
import matplotlib.colors as mcol
from matplotlib.colors import ListedColormap
import gzip

__version__ = '1.1'
__docformat__ = "reredtext en"
__what__= sys.argv[0]+": This program detects peaks in pole count maps using the Fellwalker algorithm (starlink implementation)"
#
parser = argparse.ArgumentParser(description='Plot mGC3/nGC3/GC3 pole count maps')
parser.add_argument('infile',metavar='infile',help='Input file containing pole count maps (*.cts file)',nargs=1,action='store')
parser.add_argument("-l", "--llist", action="store_true",help='Take infile as list of mgc3.cts files')
parser.add_argument('-m',help='Plot pole count map. Default is mGC3', action='store',default='mGC3',choices=['mGC3','nGC3','GC3','GC3hel','mGC3hel','smooth','usharpc','usharpn','usigma'])
parser.add_argument('-f','--fig',help='Output plot type png/eps. Default is png', action='store',default='png',choices=['png','eps','pdf'])
parser.add_argument('-proj',help='Projection npaeqd/ortho/mollweide. Default is npaeqd', action='store',default='npaeqd',choices=['npaeqd','ortho','moll'])
parser.add_argument('-log',help='Plot pole-count map in log-scale', action='store_true',default=False)
parser.add_argument('-lon0',help='Longitude for Y-axis. Default is 0.', action='store',default=0.,type=np.float)
parser.add_argument('-lat0',help='Bounding latitude for plot. Default is 90.', action='store',default=0.,type=np.float)
parser.add_argument('-dlat',help='Spacing between parallels. Default is 20.', action='store',default=20.,type=np.float)
parser.add_argument('-dlon',help='Spacing between meridians. Default is 20.', action='store',default=30.,type=np.float)
parser.add_argument('-latmax',help='Max latitude upto which meridians are drawn. Default is 80.', action='store',default=80.,type=np.float)
parser.add_argument('-mlab','--merlabels',help='Show meridian labels. Default is False', action='store_true',default=False)
parser.add_argument('-mlabr','--merlabelsr',help='Show meridian labels (right axes). Default is False', action='store_true',default=False)
parser.add_argument('-plab','--parlabels',help='Show parallel labels. Default is False', action='store_true',default=False)
parser.add_argument('-vmin',help='Min counts for color-scale. Default is min(cts)', action='store',default=None,type=np.float)
parser.add_argument('-vmax',help='Max counts for color-scale. Default is max(cts)', action='store',default=None,type=np.float)
parser.add_argument('-ms',help='Marker size. Default: 15/40 for npaeqd/ortho.', action='store',default=-1.,type=np.float)
parser.add_argument('-c','--contour',help='Plot pole-count contour map instead of raw grid.', action='store_true',default=False)
parser.add_argument('-t','--twohemispheres',help='Plot both hemispheres in pole-count map.', action='store_true',default=False)
parser.add_argument('-s','--show',help='Show plot in window. Default is False', action='store_true',default=False)
parser.add_argument('-title',help='Plot title', action='store',default=None)
parser.add_argument('-pls',metavar='PLSFILE',help='Overplot poles from peakdetect output file (.pls)', action='store',default=None)
parser.add_argument('-al','--alpha',help='Clump transparency. Default 0.4', action='store',default=0.4,type=np.float)
parser.add_argument('-ff','--ffonts',help='Increase size tick and axes labels by factor ff. Default 1.', action='store',default=1.0,type=np.float)
parser.add_argument('-flab','--flabels',help='Increase size of peak labels by factor flab. Default 1.', action='store',default=1.0,type=np.float)
parser.add_argument('-fcirc','--fcirc',help='Increase size of peak markers by factor fcirc. Default 1.', action='store',default=1.0,type=np.float)
parser.add_argument('-cmap',help='Choose color map (any matplotlib cm). Default is sron', action='store',default=None)
parser.add_argument('-ext',metavar='outfile_ext',help='Output suffix [optional]. If given output will be infile.outfile_ext.mgc3.pst',action='store',default=['',],nargs=1)

#---------Parse----------------------------
args = parser.parse_args()

#Parse inputs
if not args.llist:
 print('Reading file: %s' % (args.infile))
 file_list=[args.infile[0],]
else:
 print('Reading input files from list file: ', args.infile[0])
 file_list=scipy.genfromtxt(args.infile[0],dtype=str)
 if np.ndim(file_list)==0: file_list=np.array([file_list,])

#Mode-------------------------------------------
mode=args.m.lower()
mode_ori=args.m
print('Pole counts plotted: ', mode_ori)
if 'hel' in mode and 'mgc3' in mode:  counts_col=4-1
elif 'gc3hel' in mode: counts_col=7-1
elif 'mgc3' in mode: counts_col=3-1
elif 'ngc3' in mode: counts_col=6-1
elif 'gc3'  in mode: counts_col=5-1
elif 'smooth'  in mode: counts_col,mode_ori=3-1,'Smooth'
elif 'usharpc' in mode: counts_col,mode_ori=4-1,'Unsharp-masked'
elif 'usharpn' in mode: counts_col,mode_ori=5-1,'Unsharp-masked'
elif 'usigma'  in mode: counts_col,mode_ori=6-1,'Background St.Dev.'

#Parse raw/contour mode-------------------------
if args.contour: 
  pmode='c'
  print('Plotting contour pole-count map')
else: 
  pmode='r'
  print('Plotting raw pole-count map')
if args.log:
  print('Plotting pole-count map in log-scale')
  pmode=pmode+'l'

if args.twohemispheres:
  print('Plotting both hemispheres in pole-count map')
else: print('Plotting one hemisphere in pole-count map.')

print('Plotting using projection:', args.proj)

ori='vertical'
ori='horizontal'
ni=0

#colormap=plt.cm.jet
#if 'inferno' in args.cmap:
# colormap=newcmap.inferno
#elif 'viridis' in args.cmap:
# colormap=newcmap.viridis
#elif 'gray' in args.cmap:
# if '_r' not in args.cmap: colormap=plt.cm.gray
# else: colormap=plt.cm.gray_r
#else:
# colormap=myutils.get_sron_rainbow(N=11)

if args.cmap is None or 'sron' in args.cmap:
 colormap=myutils.get_sron_rainbow(N=11)
elif 'inferno' in args.cmap:
 colormap=newcmap.inferno
elif 'viridis' in args.cmap:
 colormap=newcmap.viridis
else:
 colormap=getattr(plt.cm,args.cmap)

if 'usharpn'  in mode:
  print('Selecting n-sigma colorlist...')
  if args.cmap is not None and 'spectral' in args.cmap:
   colormap=plt.cm.spectral
   cutcolormap=plt.cm.spectral(np.linspace(0., 1.,100)) 
   colormap = mcol.ListedColormap(cutcolormap)
   colorsel={'cmap':colormap}
  elif args.cmap is not None:
   colorsel={'cmap':colormap}
  else: #Default color list for nsigma mode
   #colorlist=['w','#FFFFE5','#FFF7BC','#b2df8a','#33a02c','#a6cee3','#1965B0','#114477','#fdbf6f','#ff7f00','#e31a1c','#771111'] # me gusta todo - verde,azul,morados,rojos
   colorlist=['w','#FFFFE5','#FFF7BC','#FEE391','#4EB265','#117733','#a6cee3','#1965B0','#114477','#ff7f00','#e31a1c','#771111']
   colorsel={'colors':colorlist}

for infilen in file_list:

  #Default title----------------------------------
  if args.title is None or args.llist: args.title=infilen

  if 'gz' in infilen: infile=gzip.open(infilen,'r')
  else: infile=open(infilen,'r')

  phio,thetao,pole_ctso=pdat=scipy.genfromtxt(infile,comments='#',usecols=(0,1,counts_col),unpack=True)

  figname_root=infilen.replace('.mgc3.cts',args.ext[0])  #works well if args.ext is empty
  figname='%s.%s.%s.%s.%s' % (figname_root,mode,args.proj[:3],pmode,args.fig)
  print('Output filename:', figname)

  phi2=np.append(phio,(phio+180.) % 360.)
  theta2=np.append(thetao,-thetao)
  pole_cts2=np.append(pole_ctso,pole_ctso) 

  if args.twohemispheres or 'moll' in args.proj:
   phi,theta,pole_cts=phi2,theta2,pole_cts2
  else: 
   phi,theta,pole_cts=phio,thetao,pole_ctso 
  #This is done this way so that for contour plotting, the twohemispheres are used for the interpolation
  #into a regular grid, but only one is shown 

  #----------------Pole count map------------------------------
  nx,ny=1,1

  mer_grid=[0.,360.,args.dlon]
  if 'npa' in args.proj: par_grid=[-args.dlat,+90.,args.dlat]
  else: par_grid=[-90.,+90.,args.dlat]

  if 'npa' in args.proj or 'moll' in args.proj:
    #For npa and moll projections, plot map as viewed from lon0 only
    if 'npa' in args.proj: fig=plt.figure(1,figsize=(8,8))
    elif 'moll' in args.proj: fig=plt.figure(1,figsize=(10,6)) 
    dw=0.8
    wo=(1.-dw)/2.
    wyo=0.05*wo
    if 'moll' in args.proj: xnudge=0.05
    else: xnudge=0.
    fig.subplots_adjust(left=wo-xnudge,right=dw+wo-xnudge,top=dw+wo,bottom=wyo)
    nrow,ncol,nplot=1,1,1
    opts=[(1,args.lon0),] 
    proj_dict={'boundinglat':args.lat0,'resolution':'l'}
    if args.ms==-1: 
      if 'npa' in args.proj: ms=30.    #tested
      elif 'moll' in args.proj: ms=12. #tested
      else: ms=70.
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
    grid_color=(0.65,0.65,0.65)
    mlabels_dic={'labels':[0,0,0,0],'labelstyle':'+/-'}
    if args.merlabels:  mlabels_dic['labels'][0]=1
    if args.merlabelsr: mlabels_dic['labels'][1]=1
    if args.parlabels: plabels_dic={'labels':[0,0,0,1],'labelstyle':'+/-'}
    else: plabels_dic={'labels':[0,0,0,0]}
    mlabels_dic['fontsize']=14.*args.ffonts
    plabels_dic['fontsize']=14.*args.ffonts
    m.drawmeridians(np.arange(mer_grid[0],mer_grid[1],mer_grid[2]),color=grid_color,linewidth=1.,
                     latmax=args.latmax,**mlabels_dic)
    m.drawparallels(np.arange(par_grid[0],par_grid[1],par_grid[2]),color=grid_color,linewidth=1.,**plabels_dic)
    m.drawmapboundary()

    x,y=m(phi,theta)

    if 'moll' in args.proj:
      #In Mollweide proj draw labels by hand
      for ptheta in np.arange(par_grid[0]+par_grid[2],par_grid[1],par_grid[2]):
        if ptheta<0.: f=0.7
        else: f=1.
        xpar,ypar=m(args.lon0-180.,ptheta)
        ax.text(f*xpar,ypar,'$%+d^\circ$' % (ptheta),horizontalalignment='right',verticalalignment='center',fontsize=10*args.ffonts)
      for mlon in np.arange(mer_grid[0],mer_grid[1],mer_grid[2])[1::2]:
        xmer,ymer=m(mlon,0.)
        ax.text(xmer,ymer,'$%d^\circ$' % (mlon),horizontalalignment='left',verticalalignment='top',fontsize=10*args.ffonts,color=grid_color)


    lmax=np.floor(np.log10(np.max(pole_cts)))
    if 'r' in pmode and 'usharpn' not in mode:
       if args.log: c=m.scatter(x,y,c=np.log10(pole_cts),edgecolor='none',s=ms,cmap=colormap,vmin=args.vmin,vmax=args.vmax)
       else:        
         if args.vmin is not None: vmin=args.vmin/10**lmax
         else: vmin=args.vmin
         if args.vmax is not None: vmax=args.vmax/10**lmax
         else: vmax=args.vmax
         c=m.scatter(x,y,c=pole_cts/10**lmax, edgecolor='none',s=ms,cmap=colormap,vmin=vmin,vmax=vmax)
    else:
       npix=500
       clevels=30
       xi = np.linspace(np.min(x),np.max(x),npix)
       yi = np.linspace(np.min(y),np.max(y),npix)
       # zi = plt.griddata(x,y,pole_cts,xi,yi, interp='linear') #,'nn') #less portable option
       grid_x, grid_y = np.meshgrid(xi,yi)
       grid_phi, grid_theta = m(grid_x, grid_y,inverse=True)
       zi = scipy.interpolate.griddata((x,y), pole_cts, (grid_x, grid_y), method='nearest')
       zi[grid_theta<0]=np.nan

       if 'usharpn'  in mode:
          sigmax=12.  #Maximum Nsigma for contour and color display
          zi_sharp_cut=zi.copy()
          zi_sharp_cut[zi_sharp_cut>=sigmax]=sigmax
          c=m.contourf(xi,yi,(zi_sharp_cut), np.arange(0.,sigmax+1.,1.),**colorsel)
       else:
        if args.log: c=m.contourf(grid_x, grid_y,np.log10(zi),clevels,vmin=args.vmin,vmax=args.vmax,**colorsel)
        else:
         if args.vmin is not None: vmin=args.vmin
         else: vmin=np.min(zi[zi>=0.]) # to avoid nan
         if args.vmax is not None: vmax=args.vmax
         else: vmax=np.max(zi[zi>=0.]) #to avoid nan
         zii=zi.copy()
         zii[(zii<vmin)]=vmin
         zii[(zii>vmax)]=vmax
         lmax=np.floor(np.log10(vmax))
         if vmax<100: lmax=0.
         c=m.contourf(grid_x, grid_y,zii/10**lmax, clevels,cmap=colormap,vmin=args.vmin,vmax=args.vmax)

    #Labels and such
    if 'npa' not in args.proj: ax.set_title('%s pole-counts' % (mode_ori))

  #If given, read in pls file
  if args.pls is not None:
   poleIDs,phis,thetas,pheight=scipy.genfromtxt(args.pls,unpack=True,usecols=(0,1,2,-1))
   #Make sure all poles are in northern hemisphere
   tmask=thetas<0.
   phis[tmask]=(phis[tmask]+180.) % 360.
   thetas[tmask]=-thetas[tmask]
   #If Molweide projection is used, plot 180deg centered on lon0
   if 'moll' in args.proj:
      phis=np.append(phis,(phis+180.) %360.)
      phis_pi=phis.copy()
      phis_pi[phis>180.]=phis_pi[phis>180.]-360.
      thetas=np.append(thetas,-thetas)
      poleIDs=np.append(poleIDs,poleIDs)
      pheight=np.append(pheight,pheight)
      if np.abs(args.lon0)<=90.: pmask=(phis_pi>=args.lon0-90.) & (phis_pi<args.lon0+90.)
      else: pmask=(phis>=(args.lon0-90.)) & (phis<(args.lon0+90.)) 
      phis,thetas,poleIDs,pheight=phis[pmask],thetas[pmask],poleIDs[pmask],pheight[pmask]
   #Keep unique pole IDs
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
     #Peak ID labels (weighed by counts)
     #u_xpeak,u_ypeak=xpoles[idmask][np.argmax(pheight[idmask])],ypoles[idmask][np.argmax(pheight[idmask])]
     u_xpeak=(xpoles[idmask]*pheight[idmask]).sum()/pheight[idmask].sum()
     u_ypeak=(ypoles[idmask]*pheight[idmask]).sum()/pheight[idmask].sum()
     m.scatter(u_xpeak,u_ypeak,c='w',alpha=0.5,edgecolor='k',s=110*(args.fcirc),zorder=100)
     ax.text(u_xpeak,u_ypeak,'%d' % (u_pid[kk]),fontsize=7*args.flabels,color='black',
                horizontalalignment='center',verticalalignment='center',zorder=101)


  #Plot colorbar
  plt.rc('xtick', labelsize=13.*args.ffonts)
  plt.rc('ytick', labelsize=13.*args.ffonts)
  if 'npa' in args.proj:
   cax0=plt.gca().get_position()
   cax=plt.axes([cax0.x0,cax0.y0+dw+0.01,dw,0.02])
   if 'usharpn' in mode:  tlocator,tformat=plt.MultipleLocator(2.),'%d'
   else: tlocator,tformat=plt.MaxNLocator(nbins=6,prune='both'),'%4.1f'
   cb=plt.colorbar(c,cax=cax,orientation='horizontal',format=tformat,ticks=tlocator)
   cax.xaxis.set_ticks_position('top')
  elif 'moll' in args.proj:
   cax0=plt.gca().get_position()
   cax=plt.axes([cax0.x0+dw+0.03,cax0.y0+(0.2),0.015,dw*0.65])
   tlocator,tformat=plt.MaxNLocator(nbins=6,prune='both'),'%4.1f'
   cb=plt.colorbar(c,cax=cax,orientation='vertical',format=tformat,ticks=tlocator)

  if 'npa' in args.proj or 'moll' in args.proj:
   #Labels and such
   if lmax>0: factorl='$\\times 10^{%d}$ ' % (lmax)
   else: factorl=''
   if 'usharpn' in mode: 
    cax.set_xlabel('%s significance ($N\sigma$)' % (mode_ori),fontsize=16.*args.ffonts,labelpad=4.*args.ffonts)
   elif args.log:
    cblabel='%s (log-stars/pole)' % (mode_ori)
    if 'npa' in args.proj: 
      cax.set_xlabel(cblabel,fontsize=16.*args.ffonts,labelpad=4.*args.ffonts)
      cax.xaxis.set_label_position('top')
    elif 'moll' in args.proj:
      cax.set_ylabel(cblabel,fontsize=16.*args.ffonts,labelpad=4.*args.ffonts) 
      cax.yaxis.set_label_position('right')
   else:
    cblabel='%s (%s stars/pole)' % (mode_ori,factorl)
    #cax.set_xlabel(cblabel,fontsize=16.*args.ffonts,labelpad=4.*args.ffonts)
    #the following lines used to be outside the else
    if 'npa' in args.proj: 
      #cax.set_xlabel(cblabel,fontsize=16.*args.ffonts,labelpad=4.*args.ffonts)
      cax.xaxis.set_label_position('top')
      cax.set_xlabel(cblabel,fontsize=16.*args.ffonts,labelpad=4.*args.ffonts)
    elif 'moll' in args.proj:
      cax.set_ylabel(cblabel,fontsize=16.*args.ffonts,labelpad=4.*args.ffonts) 
      cax.yaxis.set_label_position('right')

  if args.title:
    ax.text(0.5,1.1+0.05*args.ffonts,args.title,transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontsize=16*args.ffonts)

  fig.savefig(figname)
  if args.show: plt.show()
  else: fig.clf()


