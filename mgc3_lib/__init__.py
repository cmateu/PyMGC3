#!/usr/bin/env python
import numpy as np
import sys
import scipy 
from functools import partial
import bovy_coords as bc
import myutils
import os
import math
import gzip 
import argparse

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


'''
NAME:
   mgc3_lib
PURPOSE:
   Library containing class definitions necessary for mgc3 pole counting 
HISTORY: 
   2014-04-11 - Included optional flag to select method when extracting stars associated to poles
   2014-04-04 - Standardize for distribution
   2013-11-21 - Added ngc3 pole counts attribute to pole-grid object (mgc3, without radial velocity term)
   2013-10-xx - Fractional area per pole computation added
   2013-09-20 - Added gc3 pole counts attribute to pole-grid object
   2013-09-xx - Written. Cecilia Mateu (CIDA,IA-UNAM Ensenada)
VERSION:
   0.4 
'''

global d2r
d2r=np.pi/180.


def read_inputcat_for_mgc3(filename,pardic=None):

  #Open file
  if '.gz' in filename:
     inputfile=gzip.open(filename,'r')
     filename=filename.replace('.gz','')
  else: inputfile=open(filename,'r')

  obsdata = scipy.genfromtxt(inputfile,comments='#')

  #Do cuts
  mask = obsdata[:,0]==obsdata[:,0]  #Initialize mask to all-True-vector
  if pardic:
    for NAUX in range(1,pardic['NAUX']+1,1):
      mykey_col='AUX%d_col' % (NAUX)
      mykey_valo='AUX%d_o' % (NAUX)
      mykey_valf='AUX%d_f' % (NAUX)
      #Skip if col=998   
      if pardic[mykey_col]!=998:
       print ' Cutting input catalogue with %.1f<%s[%d]<%.1f' % (pardic[mykey_valo],mykey_col,pardic[mykey_col]+1,pardic[mykey_valf])
       #Create mask 
       mask_i = (obsdata[:,pardic[mykey_col]]>pardic[mykey_valo]) & (obsdata[:,pardic[mykey_col]]<pardic[mykey_valf])
       #Combine masks
       mask = mask & mask_i

  #Apply mask
  obsdata=obsdata[mask,:]
    
  #Return data
  return (obsdata,filename)

def equatorial2galactic_helio(alpha,delta):
   #This coordinate transformation assumes ra,dec are in the FK4 (B1950) equinox
   d2r=np.pi/180.
   lo,bo=33.*d2r,62.6*d2r
   alphao=282.25*d2r
   #Shortcuts
   sinbo,cosbo=np.sin(bo),np.cos(bo)
   sindelta,cosdelta=np.sin(delta),np.cos(delta)
   sinalpha1,cosalpha1=np.sin(alpha-alphao),np.cos(alpha-alphao)

   #Get latitude   
   sinb=sindelta*cosbo - cosdelta*sinalpha1*sinbo
   cosb=np.sqrt(1-sinb**2)  #This is ok since cosb>0 for all b in [-pi/2,+pi/2]
   b=np.arctan2(sinb,cosb)
   #Get longitute
   sinl1=(1./cosb)*(sindelta*sinbo+cosdelta*sinalpha1*cosbo)
   cosl1=(1./cosb)*(cosdelta*cosalpha1)
   #tangent of half-angle
   tanlhalf=sinl1/(1.+cosl1)
   l=2*np.arctan(tanlhalf)
   l=l % (2*np.pi)
   l=lo + l

   return (l,b)

def print_sample_parfile():
 s='''#=============================================================================================
#Paramameter file (NOTE: Column numbers are Fortran-style, i.e. 1 is the first column)
#=============================================================================================
#deg         =  If True, lat,lon in degrees. If False, radians
#lon_col     =  Longitude column
#lat_col     =  Latitude  column
#coo_glactc  =  If True, lat,lot AND corresponding proper motions 
#               assumed to be galactic (l,b), if not, assumed to be equatorial (RA,DEC)
#par_col     =  Parallax column
#par_muas    =  If True units for parallax assumed to be muas, if False mas
#---------------------------------
#pm_lon      =  Column for proper motion in longitude direction
#pm_lon_red  =  If True, mu_l is reduced proper motion (i.e. mu_l*cos(b))
#pm_lat_col  =  Column for proper motion in latitude direction
#vrad_col    =  Radial Velocity column
#pm_muas     =  If True units for proper motions assumed to be muas/yr, if False mas/yr
#---------------------------------
#tol_r       =  r-tolerance for mgc3 pole-counting
#tol_v       =  v-tolerance for mgc3 pole-counting
#tol_deg     =  If True, tolerances are assumed to be in degrees. If False, radians
#grid_step   =  Step (in same units as tol) for pole grid
#grid_lon_o  =  Initial lon for pole grid
#grid_lon_f  =  Final lon for pole grid
#grid_lat_o  =  Initial lat for pole grid
#grid_lat_f  =  Final lat for pole grid
#---------------------------------
#AUX1_col    =  Auxiliary column. Catalogue ata with AUX1_o<AUX1<AUX1_f will be selected
#AUX1_o      =  Any number of AUX? cols can be used.
#AUX1_f      =  For multiple AUX? columns, the criteria is combined with boolean AND 
#----------------Coordinate params------------------
deg         False
lon_col     1
lat_col     2
coo_glactc  True
par_col     3
par_muas    True
#----------------Proper motion params---------------
pm_lon_col  4
pm_lon_red  True
pm_lat_col  5
pm_muas     True
vrad_col    6
#----------------mGC3 params------------------------
tol_r       2.
tol_v       2.
tol_deg     True
grid_step   2.0   #All grid parameters should be in degrees
grid_lon_o  0.
grid_lon_f  360.
grid_lat_o  0.
grid_lat_f  90.
#---------------Additional pars---------------------
#AUX1_col    7    #Auxiliary column. Only catalogue data with AUX1_o<AUX1_col<AUX1_f will be used
#AUX1_o      0.   #Any number of AUX? cols can be used.
#AUX1_f      20.  #For multiple AUX? columns, the criteria is combined with boolean AND 
#----------------end------------------------'''

 ofile=open('mgc3_sample.par','w')
 ofile.write(s)
 ofile.close()

class print_parfile_action(argparse.Action):
     def __call__(self, parser, namespace, values, option_string=None):
            print_sample_parfile()
            sys.exit('Sample parfile mgc3_sample.par created')

#Parse parameter file
def parse_pars(parfile):
  pf=scipy.genfromtxt(parfile,comments='#',dtype='S')

  naux,dic=0,{}
  for line in pf:
     if 'col' in line[0]: 
        dic[line[0]]=int(line[1])-1 #So 1st column corresponds to 0 (input is human or F-style)
        if 'AUX' in line[0]: naux=naux+1  #Count how many auxiliary criteria are defined
     elif 'tol'   in line[0] and 'deg' not in line[0]: dic[line[0]]=float(line[1])
     elif 'grid_' in line[0]: dic[line[0]]=float(line[1])
     elif 'AUX'   in line[0]: dic[line[0]]=float(line[1])
     else: dic[line[0]]= line[1].lower() == 'true'  #Convert to boolean
  
  #Save number of auxialiry cols in dictionary
  dic['NAUX']=naux

  return dic

#Global constants 
def load_constants():
  global rsun
  global Ugc_hel
  global Vgc_hel
  global Wgc_hel
  global vsun
  global d2r
  global Ap
  global Av
  rsun=8.5          #In kpc
  Ugc_hel=10.3      #In km/s
  Vgc_hel=232.6     #In km/s
  Wgc_hel=5.9       #In km/s
  vsun=np.sqrt(Ugc_hel**2 + Vgc_hel**2 + Wgc_hel**2)
  d2r=np.pi/180.
  Ap=1000.    #muas.kpc
  Av=4.74047 #yr km/s

class my_constants(object):
                         #kpc         km/s         km/s         km/s  muas.kpc   yr km/s
   def __init__(self,rsun=8.5,Ugc_hel=10.3,Vgc_hel=232.6,Wgc_hel=5.9,Ap=1000.,Av=4.74047):

     self.rsun=rsun
     self.Ugc_hel=Ugc_hel
     self.Vgc_hel=Vgc_hel
     self.Wgc_hel=Wgc_hel
     self.vsun=np.sqrt(self.Ugc_hel**2 + self.Vgc_hel**2 + self.Wgc_hel**2)
     self.d2r=np.pi/180.
     self.Ap=Ap
     self.Av=Av

#---------------------------------------------------------------------------------------
# ini_polegrid: Iinitialize pole grid
#---------------------------------------------------------------------------------------

class pole_grid(my_constants):

  '''This class creates a pole-grid object with pole-counts initialized to zero

    Parameters
    ----------
    poles : float or [array, array], optional
      Pole grid parameters:

      * If float, it is the grid spacing in degrees 
      * If [array,array], pole grid longitudes and latitudes in degrees (longs,lats=poles)

    cst : my_constants object instance - optional. my_constants attributes are inherited by this class.
          This objects contains the values for relevant constants for
          mGC3 calculations, i.e. rsun, (U,V,W)_GC, Ap, Av. If not explicitly provided by the user
          it will take an instance of the my_constance class using the default values

    Returns
    -------

    P : object of pole-grid class
'''

  def __init__(self,poles=2.,cst=None,pole_grid_dic=None):

   #Inherit my_constance object attributes 
   if cst:
     my_constants.__init__(self,rsun=cst.rsun,Ugc_hel=cst.Ugc_hel,Vgc_hel=cst.Vgc_hel,Wgc_hel=Wgc_hel,Ap=cst.Ap,Av=cst.Av)
   else:
     my_constants.__init__(self)

   #Initialize empty arrays
   self.l,self.b=np.array([]),np.array([])
   self.mgc3hel,self.np_mgc3, self.np_gc3=np.array([]),np.array([]),np.array([])
   self.np_ngc3,self.farea=np.array([]),np.array([])
   self.np_mgc3_mask,self.np_gc3_mask=np.array([]),np.array([])
   self.sinlp,self.coslp,self.sinbp,self.cosbp=np.array([]),np.array([]),np.array([]),np.array([])
   self.L_rsun,self.L_vsun=np.array([]),np.array([])

   #Cycle through grid
   if np.ndim(poles)==0:
     if pole_grid_dic: 
        lo,lf=pole_grid_dic['grid_lon_o'],pole_grid_dic['grid_lon_f']
        bo,bf=pole_grid_dic['grid_lat_o'],pole_grid_dic['grid_lat_f']
     else: lo,lf,bo,bf=0.,360.,0.,90.
     gstep=poles
     print 'Building grid with gstep=%.2f, (%.2f<lon<%.2f,%.2f<lat<%.2f)' % (gstep,lo,lf,bo,bf)
     for b in np.arange(bo,bf,gstep):
       lstep=gstep/np.cos(b*d2r)
       for l in np.arange(lo,lf,lstep):
         self.ini_pole_count_single(l,b)       
   elif np.ndim(poles)==1:
      ls,bs=poles
      print 'Using single pole:', ls,bs
      self.ini_pole_count_single(ls,bs)
   elif np.ndim(poles)==2:
      ls,bs=poles
      print 'Building grid with arbitrary poles:', ls,bs
      for l,b in ls,bs:
         self.ini_pole_count_single(l,b)
   else:
     sys.exit('Wrong dimension of object passed to poles keyword in pole_grid class. Input dim is: %d' % (np.ndim(poles)))
    
  def ini_pole_count_single(self,l_deg,b_deg):
     lr,br=l_deg*d2r,b_deg*d2r
     self.l,self.b=np.append(self.l,l_deg),np.append(self.b,b_deg)
     #Count attributes
     self.mgc3hel=np.append(self.mgc3hel,0.)
     self.np_mgc3=np.append(self.np_mgc3,0.)
     self.np_gc3=np.append(self.np_gc3,0.)
     self.np_ngc3=np.append(self.np_ngc3,0.)
     self.farea=np.append(self.farea,0.)
     #Mask attributes
     self.np_mgc3_mask=np.append(self.np_mgc3_mask,False)
     self.np_gc3_mask=np.append(self.np_gc3_mask,False)
     self.sinlp,self.coslp=np.append(self.sinlp,np.sin(lr)),np.append(self.coslp,np.cos(lr))
     self.sinbp,self.cosbp=np.append(self.sinbp,np.sin(br)),np.append(self.cosbp,np.cos(br))
     self.L_rsun=np.append(self.L_rsun, -self.rsun*np.cos(lr)*np.cos(br))
     self.L_vsun=np.append(self.L_vsun, self.Ugc_hel*np.cos(lr)*np.cos(br) + self.Vgc_hel*np.sin(lr)*np.cos(br) + self.Wgc_hel*np.sin(br))

  #---------------------------------------------------------------------------------------
  # modGC3 pole count method
  #---------------------------------------------------------------------------------------
  def mgc3(self, obsdata, pars=None, return_mask=False):

   #Pass fixed keyword values
   map_modgc3 = partial(self.mgc3_single_obs, pars=pars, return_mask=return_mask)
   #Loop over observations
   map(map_modgc3, obsdata)

   #Normalize appropriately auxiliary counts for area fraction computation
   if pars['tol_deg']: tolr=pars['tol_deg']*self.d2r,
   fgc=np.sin(tolr)  #The fractional area of the sphere covered by any great circle is 4pi*sin(tolr)/4pi
   self.farea=self.farea/(fgc*obsdata[:,0].size)

  def mgc3_allobs_one_pole(self, obsdata, pars=None, return_mask=False):

   #This function call will act on the full observations object
   #NOTE: obsdata needs to be transposed so obsdata[col_n] will pick out the n-th column
   #rather than the n-row (which would be default behaviour, since for map-compatibility,
   #we're not using obs[:,col_n] but rather obs[col_n])
   pole_mask=self.mgc3_single_obs(obsdata.T,pars=pars,return_mask=return_mask)

   return pole_mask

  def mgc3_single_obs(self, obs, pars=None, return_mask=False):
  
    '''This method does mGC3 pole-counts (Mateu et al. 2011) in a given pole-grid, for the given observational data
  
      Inputs
      ----------
      obs :   2D-array containing observational data (i.e. l,b,parallax,vrad,pm_lon,pm_lat)
  
      Parameters
      ----------
      pars  : Parameters dictionary - optional. This defines which columns contain each of the observables.
  	    Default value is None. In this case, a default parameter file is printed and default
  	    values are assumed for each parameter. 
  
      Returns
      ----------
      pgrid : Updates the input pole_grid object counts (attributes pgrid.np_mgc3 and pgrid.np_gc3)
    '''

    if pars is None:
     print 'No params file found. Using default parameter file mgc3.par'
     print_sample_parfile()
     pars=parse_pars('mgc3.par')

    tolr=pars['tol_r']
    tolv=pars['tol_v']
    if pars['tol_deg']: tolr,tolv=tolr*d2r,tolv*d2r
    sin_tolr=np.sin(tolr)
    sin_tolv=np.sin(tolv)
  
    #Check coord system. If equatorial, convert coords and proper motions to galactic
    #Makes use of Jo Bovy's coordinate conversion library
    if pars['deg']: fc=d2r 
    else: fc=1.
    #fc will be used to convert inputs to radians from the start (if necessary)
    if not pars['coo_glactc']:
      ra,dec=fc*obs[pars['lon_col']],fc*obs[pars['lat_col']]
      pmra_red,pmdec=obs[pars['pm_lon_col']],obs[pars['pm_lat_col']]
      #Reduce pm_lon if necessary 
      if not pars['pm_lon_red']: pmra_red=pmra_red*np.cos(dec)
      #Convert. Output is in the same units as input
      if np.ndim(ra)>0:  
        lb_m=bc.radec_to_lb(ra,dec,degree=False,epoch=2000.0)
        l,b=np.array(list(lb_m[:,0])),np.array(list(lb_m[:,1])) #Yes, this is awful,but if not these are not proper arrays
        mm=bc.pmrapmdec_to_pmllpmbb(pmra_red,pmdec,ra,dec,degree=False,epoch=2000.0)
        pml_red,pmb=np.array(list(mm[:,0])),np.array(list(mm[:,1]))
      else: 
        l,b=bc.radec_to_lb(ra,dec,degree=False,epoch=2000.0) 
        pml_red,pmb=bc.pmrapmdec_to_pmllpmbb(pmra_red,pmdec,ra,dec,degree=False,epoch=2000.0)
    else:
      l,b=fc*obs[pars['lon_col']],fc*obs[pars['lat_col']]
      pml_red,pmb=obs[pars['pm_lon_col']],obs[pars['pm_lat_col']]
      if not pars['pm_lon_red']: pml_red=pml_red*np.cos(b)
  
    #Compute trig functions
    cos_b,sin_b=np.cos(b),np.sin(b)  
    cos_l,sin_l=np.cos(l),np.sin(l)  
  
    #Read parallax and radial velocity
    if pars['par_muas']: factor=1.
    else: factor=1000.
    parallax=factor*obs[pars['par_col']]
    vrad=obs[pars['vrad_col']]
  
    #Angular velocities
    if pars['pm_muas']: factor=1.
    else: factor=1000.
    omega_r=parallax*vrad                 #parallax*Vr
    omega_l=self.Av*pml_red*factor      #4.74*mu_l* <- must use mul*=mul*cos(b). This has been triple-checked (last triple check 01/2014)
    omega_b=self.Av*pmb*factor      #4.74*mu_b
  
    #Compute vhel, rhel
    vhel=np.sqrt(omega_r**2 + omega_l**2 + omega_b**2)
    rhel=self.Ap/parallax
  
    #Trig-functions for pole grid points
    sin_lp,cos_lp,sin_bp,cos_bp=self.sinlp,self.coslp,self.sinbp,self.cosbp  

    #Producto escalar con vector posicion heliocentrico (helpos) 
    helpos_dotprod=cos_b*cos_bp*(cos_l*cos_lp + sin_l*sin_lp)+sin_b*sin_bp
    vxhel=(omega_r*cos_l*cos_b - omega_l*sin_l - omega_b*cos_l*sin_b)
    vyhel=(omega_r*sin_l*cos_b + omega_l*cos_l - omega_b*sin_l*sin_b)
    vzhel=(omega_r*sin_b + omega_b*cos_b)
    helvel_dotprod=(cos_lp*cos_bp*vxhel + sin_lp*cos_bp*vyhel + sin_bp*vzhel)

    #Producto escalar, omitiendo el termino de velocidad radial
    muxhel=-omega_l*sin_l - omega_b*cos_l*sin_b
    muyhel=omega_l*cos_l - omega_b*sin_l*sin_b
    muzhel=omega_b*cos_b
    helmu_dotprod=(cos_lp*cos_bp*muxhel + sin_lp*cos_bp*muyhel + sin_bp*muzhel)
   
    #Criterios de Posicion y Velocidad Angular Heliocentrica
    mask_poshel=np.abs(helpos_dotprod)<=sin_tolr
    mask_velhel=np.abs(helvel_dotprod/vhel)<=sin_tolv
    mask_hel=mask_poshel & mask_velhel

    if len(sin_lp)==1:
      self.mgc3hel=np.sum(1*mask_hel)
      #For f_area computation only
      self.farea=np.sum(1*mask_poshel)
    else:
      self.mgc3hel[mask_hel]=self.mgc3hel[mask_hel]+1
      #For f_area computation only
      self.farea[mask_poshel]=self.farea[mask_poshel]+1
  
    #Galactocentrica
    pirgal=np.sqrt(self.Ap**2 + (parallax*self.rsun)**2 - 2*self.rsun*self.Ap*parallax*cos_b*cos_l)
    galpos_dotprod=(self.Ap*helpos_dotprod + parallax*self.L_rsun)
    omegagal=np.sqrt((parallax*self.Ugc_hel + vxhel)**2+(parallax*self.Vgc_hel + vyhel)**2+(parallax*self.Wgc_hel + vzhel)**2)
    galvel_dotprod=(helvel_dotprod + parallax*self.L_vsun)
  
    #Galactocentrica sin velocidad radial
    mugal=np.sqrt((parallax*self.Ugc_hel + muxhel)**2+(parallax*self.Vgc_hel + muyhel)**2+(parallax*self.Wgc_hel + muzhel)**2)
    galmu_dotprod=(helmu_dotprod + parallax*self.L_vsun)

    #Criterios de Posicion y Velocidad Angular Galactocentrica
    mask_posgal=np.abs(galpos_dotprod)<=sin_tolr*pirgal
    mask_velgal=np.abs(galvel_dotprod)<=sin_tolv*omegagal
    mask_gal=mask_posgal & mask_velgal
    #Criterio sin velocidad radial
    mask_mugal=np.abs(galmu_dotprod)<=sin_tolv*mugal
    mask_muposgal=mask_posgal & mask_mugal

    if len(sin_lp)==1:
      self.np_mgc3=np.sum(1*mask_gal)
      self.np_gc3=np.sum(1*mask_posgal)
      self.np_ngc3=np.sum(1*mask_muposgal)
    else:
      self.np_mgc3[mask_gal]=self.np_mgc3[mask_gal]+1
      self.np_gc3[mask_posgal]=self.np_gc3[mask_posgal]+1 
      self.np_ngc3[mask_muposgal]=self.np_ngc3[mask_muposgal]+1

    if return_mask:
      if 'mGC3' in return_mask: 
        print '   Selecting stars fulfilling mGC3 criteria'
        return mask_gal
      elif 'nGC3' in return_mask: 
        print '   Selecting stars fulfilling nGC3 criteria'
        return mask_muposgal
      else: 
        print '   Selecting stars fulfilling GC3 criteria'
        return mask_posgal

  def get_phi_theta_for_survey(self,obs,pars=None):

    if pars is None:
     print 'No params file found. Using default parameter file mgc3.par'     
     print_sample_parfile()
     pars=parse_pars('mgc3.par')

    #Check coord system. If equatorial, convert coords and proper motions to galactic
    #Makes use of Jo Bovy's coordinate conversion library
    if pars['deg']: fc=d2r
    else: fc=1.
    #fc will be used to convert inputs to radians from the start (if necessary)
    if not pars['coo_glactc']:
      ra,dec=fc*obs[:,pars['lon_col']],fc*obs[:,pars['lat_col']]
      #Convert. Output is in the same units as input
      if np.ndim(ra)>0: 
        lb_m=bc.radec_to_lb(ra,dec,degree=False,epoch=2000.0)
        l,b=np.array(list(lb_m[:,0])),np.array(list(lb_m[:,1]))
      else: l,b=bc.radec_to_lb(ra,dec,degree=False,epoch=2000.0)
    else:
      l,b=fc*obs[:,pars['lon_col']],fc*obs[:,pars['lat_col']]

    #Read parallax and radial velocity
    if pars['par_muas']: factor=1.
    else: factor=1000.
    parallax=factor*obs[:,pars['par_col']]
    rhel=self.Ap/parallax

    xyz_hel=bc.lbd_to_XYZ(l,b,rhel,degree=False)
    xhel,ygal,zgal=xyz_hel[:,0],xyz_hel[:,1],xyz_hel[:,2]
    print np.shape(xhel)
    xgal=xhel-self.rsun

    #convert to spherical galactocentric
    Rproy=np.sqrt(xgal**2+ygal**2)
    sinphi,cosphi=ygal/Rproy,xgal/Rproy
    tanphihalf=sinphi/(1.+cosphi)
    phi=2*np.arctan(tanphihalf)
    phi=phi % (2*np.pi)
    theta=np.arcsin(zgal/Rproy)
    Rgal=np.sqrt(Rproy**2+zgal**2)

    #If deg_flag is set, convert phi,theta to radians to be consistent
    if pars['deg']: phi,theta=phi/self.d2r,theta/self.d2r

    #Add this columns to the obsdata matrix and dictionary
    new_cols=np.array([phi,theta,Rgal]).T  #matrix with p,t,r as columns

    print np.shape(obs), np.shape(new_cols)
    new_obs=np.hstack((obs,new_cols))      #append as new columns
    print np.shape(new_obs)
    Nobs=len(obs[0,:])
    new_pars=pars
    new_pars['phi_col']=Nobs
    new_pars['theta_col']=Nobs+1
    new_pars['Rgal_col']=Nobs+2

    return (new_obs,new_pars)

  def get_uniform_survey_footprint(self,obs,pars=None,c='k',ms=3.,show=True,Npts=1e4):
    if pars is None:
     print 'No params file found. Using default parameter file mgc3.par'
     print_sample_parfile()
     pars=parse_pars('mgc3.par')

    #Compute phi and theta for real survey data
    new_obs,new_pars=self.get_phi_theta_for_survey(obs,pars=pars)
    obs_lon,obs_lat=new_obs[:,pars['phi_col']],new_obs[:,pars['theta_col']]
    scipy.savetxt('out.phitheta.dat',new_obs)

    #Uniform survey
    lon,lat=myutils.uniform_spherical_data(Npts)
    lon,lat=lon/self.d2r,lat/self.d2r
    #Match to observations
    if pars['deg']: f=1.
    else: f=1./self.d2r

    scipy.savetxt('.aux1',np.array([lon,lat]).T)
    scipy.savetxt('.aux2',np.array([f*obs_lon,f*obs_lat]).T)
    
    #Match
    stilts_path='/Applications/TOPCAT.app/'
    os.system('%s/stilts -Xmx1024M -disk tskymatch2 in1=.aux1 in2=.aux2 ifmt1=ascii ifmt2=ascii ofmt=ascii out=.match ra1=col1 dec1=col2 ra2=col1 dec2=col2 error=1800 join=1and2 find=best' % (stilts_path))

    foot_lon,foot_lat=scipy.genfromtxt('.match',unpack=True,usecols=(0,1))
    print 'Matches',foot_lon.size

    #Create new survey object
    dummy=np.zeros_like(foot_lon)
    foot_survey=np.array([foot_lon,foot_lat,dummy]).T
    foot_pars=pars
    #This synthetic uniform survey is only used to compute F_area for each pole. Only position angles are necessary  
    #thats why parallax,proper motions and radial velocity are the same dummy column
    foot_pars['lon_col']=0
    foot_pars['lat_col']=1
    foot_pars['deg']=True
    foot_pars['coo_glactc']=True
    foot_pars['par_col']=2
    foot_pars['pm_lon_col']=2
    foot_pars['pm_lat_col']=2
    foot_pars['vrad_col']=2

    return (foot_survey,foot_pars)

#---------------------------------------------------------------------------------------
