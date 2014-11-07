###############################################################################
#
#   bovy_coords: module for coordinate transformations between the equatorial
#                and Galactic coordinate frame
#
#
#      Main included functions:
#            radec_to_lb
#            lb_to_radec
#            lbd_to_XYZ
#            XYZ_to_lbd
#            rectgal_to_sphergal
#            sphergal_to_rectgal
#            vrpmllpmbb_to_vxvyvz
#            vxvyvz_to_vrpmllpmbb
#            pmrapmdec_to_pmllpmbb
#            pmllpmbb_to_pmrapmdec
#            cov_pmrapmdec_to_pmllpmbb
#            cov_dvrpmllbb_to_vxyz
#            XYZ_to_galcenrect
#            XYZ_to_galcencyl
#            galcenrect_to_XYZ
#            galcencyl_to_XYZ
#            rect_to_cyl
#            cyl_to_rect
#            rect_to_cyl_vec
#            cyl_to_rect_vec
#            vxvyvz_to_galcenrect
#            vxvyvz_to_galcencyl
#            galcenrect_to_vxvyvz
#            galcencyl_to_vxvyvz
#            dl_to_rphi_2d
#            rphi_to_dl_2d
#            Rz_to_coshucosv
#            Rz_to_uv
#            uv_to_Rz
#
##############################################################################
#############################################################################
#Copyright (c) 2010 - 2012, Jo Bovy
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
from functools import wraps
import math as m
import numpy as nu
import scipy as sc
_DEGTORAD= m.pi/180.
_K=4.74047
def scalarDecorator(func):
    """Decorator to return scalar outputs as a set"""
    @wraps(func)
    def scalar_wrapper(*args,**kwargs):
        if nu.array(args[0]).shape == ():
            scalarOut= True
            newargs= ()
            for ii in range(len(args)):
                newargs= newargs+(nu.array([args[ii]]),)
            args= newargs
        else:
            scalarOut= False
        result= func(*args,**kwargs)
        if scalarOut:
            out= ()
            for ii in range(result.shape[1]):
                out= out+(result[0,ii],)
            return out
        else:
            return result
    return scalar_wrapper
def degreeDecorator(inDegrees,outDegrees):
    """Decorator to transform angles from and to degrees if necessary"""
    def wrapper(func):
        @wraps(func)
        def wrapped(*args,**kwargs):
            if kwargs.get('degree',False):
                newargs= ()
                for ii in range(len(args)):
                    if ii in inDegrees:
                        newargs= newargs+(args[ii]*nu.pi/180.,)
                    else:
                        newargs= newargs+(args[ii],)
                args= newargs
            out= func(*args,**kwargs)
            if kwargs.get('degree',False):
                for indx in outDegrees:
                    out[:,indx]/= nu.pi/180.
            return out
        return wrapped
    return wrapper            

@scalarDecorator
@degreeDecorator([0,1],[0,1])
def radec_to_lb(ra,dec,degree=False,epoch=2000.0):
    """
    NAME:

       radec_to_lb

    PURPOSE:

       transform from equatorial coordinates to Galactic coordinates

    INPUT:

       ra - right ascension

       dec - declination

       degree - (Bool) if True, ra and dec are given in degree and l and b will be as well

       epoch - epoch of ra,dec (right now only 2000.0 and 1950.0 are supported)

    OUTPUT:

       l,b

       For vector inputs [:,2]

    HISTORY:

       2009-11-12 - Written - Bovy (NYU)

       2014-06-14 - Re-written w/ numpy functions for speed and w/ decorators for beauty - Bovy (IAS)

    """
    #First calculate the transformation matrix T
    theta,dec_ngp,ra_ngp= get_epoch_angles(epoch)
    T= sc.dot(sc.array([[sc.cos(theta),sc.sin(theta),0.],[sc.sin(theta),-sc.cos(theta),0.],[0.,0.,1.]]),sc.dot(sc.array([[-sc.sin(dec_ngp),0.,sc.cos(dec_ngp)],[0.,1.,0.],[sc.cos(dec_ngp),0.,sc.sin(dec_ngp)]]),sc.array([[sc.cos(ra_ngp),sc.sin(ra_ngp),0.],[-sc.sin(ra_ngp),sc.cos(ra_ngp),0.],[0.,0.,1.]])))
    #Whether to use degrees and scalar input is handled by decorators
    XYZ= nu.array([nu.cos(dec)*nu.cos(ra),
                   nu.cos(dec)*nu.sin(ra),
                   nu.sin(dec)])
    galXYZ= nu.dot(T,XYZ)
    b= nu.arcsin(galXYZ[2])
    l= nu.arctan2(galXYZ[1]/sc.cos(b),galXYZ[0]/sc.cos(b))
    l[l<0.]+= 2.*nu.pi
    out= nu.array([l,b])
    return out.T

@scalarDecorator
@degreeDecorator([0,1],[0,1])
def lb_to_radec(l,b,degree=False,epoch=2000.0):
    """
    NAME:

       lb_to_radec

    PURPOSE:

       transform from Galactic coordinates to equatorial coordinates

    INPUT:

       l - Galactic longitude

       b - Galactic lattitude

       degree - (Bool) if True, l and b are given in degree and ra and dec will be as well

       epoch - epoch of target ra,dec (right now only 2000.0 and 1950.0 are supported)

    OUTPUT:

       ra,dec

       For vector inputs [:,2]

    HISTORY:

       2010-04-07 - Written - Bovy (NYU)

       2014-06-14 - Re-written w/ numpy functions for speed and w/ decorators for beauty - Bovy (IAS)

    """
    #First calculate the transformation matrix T'
    theta,dec_ngp,ra_ngp= get_epoch_angles(epoch)
    T= sc.dot(sc.array([[sc.cos(ra_ngp),-sc.sin(ra_ngp),0.],[sc.sin(ra_ngp),sc.cos(ra_ngp),0.],[0.,0.,1.]]),sc.dot(sc.array([[-sc.sin(dec_ngp),0.,sc.cos(dec_ngp)],[0.,1.,0.],[sc.cos(dec_ngp),0.,sc.sin(dec_ngp)]]),sc.array([[sc.cos(theta),sc.sin(theta),0.],[sc.sin(theta),-sc.cos(theta),0.],[0.,0.,1.]])))
    #Whether to use degrees and scalar input is handled by decorators
    XYZ= nu.array([nu.cos(b)*nu.cos(l),
                   nu.cos(b)*nu.sin(l),
                   nu.sin(b)])
    eqXYZ= nu.dot(T,XYZ)
    dec= nu.arcsin(eqXYZ[2])
    ra= nu.arctan2(eqXYZ[1],eqXYZ[0])
    ra[ra<0.]+= 2.*nu.pi
    return nu.array([ra,dec]).T

@scalarDecorator
@degreeDecorator([0,1],[])
def lbd_to_XYZ(l,b,d,degree=False):
    """
    NAME:

       lbd_to_XYZ

    PURPOSE:

       transform from spherical Galactic coordinates to rectangular Galactic coordinates (works with vector inputs)

    INPUT:

       l - Galactic longitude (rad)

       b - Galactic lattitude (rad)

       d - distance (arbitrary units)

       degree - (bool) if True, l and b are in degrees

    OUTPUT:

       [X,Y,Z] in whatever units d was in

       For vector inputs [:,3]

    HISTORY:

       2009-10-24- Written - Bovy (NYU)

       2014-06-14 - Re-written w/ numpy functions for speed and w/ decorators for beauty - Bovy (IAS)

    """
    #Whether to use degrees and scalar input is handled by decorators
    return nu.array([d*nu.cos(b)*nu.cos(l),
                     d*nu.cos(b)*nu.sin(l),
                     d*nu.sin(b)]).T

def rectgal_to_sphergal(X,Y,Z,vx,vy,vz,degree=False):
    """
    NAME:

       rectgal_to_sphergal

    PURPOSE:

       transform phase-space coordinates in rectangular Galactic coordinates to spherical Galactic coordinates (can take vector inputs)

    INPUT:

       X - component towards the Galactic Center (kpc)

       Y - component in the direction of Galactic rotation (kpc)

       Z - component towards the North Galactic Pole (kpc)

       vx - velocity towards the Galactic Center (km/s)

       vy - velocity in the direction of Galactic rotation (km/s)

       vz - velocity towards the North Galactic Pole (km/s)

       degree - (Bool) if True, return l and b in degrees

    OUTPUT:

       (l,b,d,vr,pmll,pmbb) in (rad,rad,kpc,km/s,mas/yr,mas/yr)

    HISTORY:

       2009-10-25 - Written - Bovy (NYU)

    """
    lbd= XYZ_to_lbd(X,Y,Z,degree=degree)
    vrpmllpmbb= vxvyvz_to_vrpmllpmbb(vx,vy,vz,X,Y,Z,XYZ=True)
    if sc.array(X).shape == ():
        return sc.array([lbd[0],lbd[1],lbd[2],vrpmllpmbb[0],vrpmllpmbb[1],vrpmllpmbb[2]])
    else:
        out=sc.zeros((len(X),6))
        out[:,0:3]= lbd
        out[:,3:6]= vrpmllpmbb
        return out

def sphergal_to_rectgal(l,b,d,vr,pmll,pmbb,degree=False):
    """
    NAME:

       sphergal_to_rectgal

    PURPOSE:

       transform phase-space coordinates in spherical Galactic coordinates to rectangular Galactic coordinates (can take vector inputs)

    INPUT:

       l - Galactic longitude (rad)

       b - Galactic lattitude (rad)

       d - distance (kpc)

       vr - line-of-sight velocity (km/s)

       pmll - proper motion in the Galactic longitude direction (mu_l*cos(b) ) (mas/yr)

       pmbb - proper motion in the Galactic lattitude (mas/yr)

       degree - (bool) if True, l and b are in degrees

    OUTPUT:

       (X,Y,Z,vx,vy,vz) in (kpc,kpc,kpc,km/s,km/s,km/s)

    HISTORY:

       2009-10-25 - Written - Bovy (NYU)

    """
    XYZ= lbd_to_XYZ(l,b,d,degree=degree)
    vxvyvz= vrpmllpmbb_to_vxvyvz(vr,pmll,pmbb,l,b,d,XYZ=False,degree=degree)
    if sc.array(l).shape == ():
        return sc.array([XYZ[0],XYZ[1],XYZ[2],vxvyvz[0],vxvyvz[1],vxvyvz[2]])
    else:
        out=sc.zeros((len(l),6))
        out[:,0:3]= XYZ
        out[:,3:6]= vxvyvz
        return out

@scalarDecorator
@degreeDecorator([3,4],[])
def vrpmllpmbb_to_vxvyvz(vr,pmll,pmbb,l,b,d,XYZ=False,degree=False):
    """
    NAME:

       vrpmllpmbb_to_vxvyvz

    PURPOSE:

       Transform velocities in the spherical Galactic coordinate frame to the rectangular Galactic coordinate frame (can take vector inputs)

    INPUT:

       vr - line-of-sight velocity (km/s)

       pmll - proper motion in the Galactic longitude (mu_l * cos(b))(mas/yr)

       pmbb - proper motion in the Galactic lattitude (mas/yr)

       l - Galactic longitude

       b - Galactic lattitude

       d - distance (kpc)

       XYZ - (bool) If True, then l,b,d is actually X,Y,Z (rectangular Galactic coordinates)

       degree - (bool) if True, l and b are in degrees

    OUTPUT:

       (vx,vy,vz) in (km/s,km/s,km/s)

       For vector inputs [:,3]

    HISTORY:

       2009-10-24 - Written - Bovy (NYU)

       2014-06-14 - Re-written w/ numpy functions for speed and w/ decorators for beauty - Bovy (IAS)

    """
    #Whether to use degrees and scalar input is handled by decorators
    if XYZ: #undo the incorrect conversion that the decorator did
        if degree:
            l*= 180./nu.pi 
            b*= 180./nu.pi 
        lbd= XYZ_to_lbd(l,b,d,degree=False)
        l= lbd[:,0]
        b= lbd[:,1]
        d= lbd[:,2]
    R=nu.zeros((3,3,len(l)))
    R[0,0]= nu.cos(l)*nu.cos(b)
    R[1,0]= -nu.sin(l)
    R[2,0]= -nu.cos(l)*nu.sin(b)
    R[0,1]= nu.sin(l)*nu.cos(b)
    R[1,1]= nu.cos(l)
    R[2,1]= -nu.sin(l)*nu.sin(b)
    R[0,2]= nu.sin(b)
    R[2,2]= nu.cos(b)
    invr= nu.array([[vr,vr,vr],
                    [d*pmll*_K,d*pmll*_K,d*pmll*_K],
                    [d*pmbb*_K,d*pmbb*_K,d*pmbb*_K]])
    return (R.T*invr.T).sum(-1)

@scalarDecorator
@degreeDecorator([3,4],[])
def vxvyvz_to_vrpmllpmbb(vx,vy,vz,l,b,d,XYZ=False,degree=False):
    """
    NAME:

       vxvyvz_to_vrpmllpmbb

    PURPOSE:

       Transform velocities in the rectangular Galactic coordinate frame to the spherical Galactic coordinate frame (can take vector inputs)

    INPUT:

       vx - velocity towards the Galactic Center (km/s)

       vy - velocity in the direction of Galactic rotation (km/s)

       vz - velocity towards the North Galactic Pole (km/s)

       l - Galactic longitude

       b - Galactic lattitude

       d - distance (kpc)

       XYZ - (bool) If True, then l,b,d is actually X,Y,Z (rectangular Galactic coordinates)

       degree - (bool) if True, l and b are in degrees

    OUTPUT:

       (vr,pmll,pmbb) in (km/s,mas/yr,mas/yr); pmll = mu_l * cos(b)

       For vector inputs [:,3]

    HISTORY:

       2009-10-24 - Written - Bovy (NYU)

       2014-06-14 - Re-written w/ numpy functions for speed and w/ decorators for beauty - Bovy (IAS)

    """
    #Whether to use degrees and scalar input is handled by decorators
    if XYZ: #undo the incorrect conversion that the decorator did
        if degree:
            l*= 180./nu.pi 
            b*= 180./nu.pi 
        lbd= XYZ_to_lbd(l,b,d,degree=False)
        l= lbd[:,0]
        b= lbd[:,1]
        d= lbd[:,2]
    R=nu.zeros((3,3,len(l)))
    R[0,0]= nu.cos(l)*nu.cos(b)
    R[0,1]= -nu.sin(l)
    R[0,2]= -nu.cos(l)*nu.sin(b)
    R[1,0]= nu.sin(l)*nu.cos(b)
    R[1,1]= nu.cos(l)
    R[1,2]= -nu.sin(l)*nu.sin(b)
    R[2,0]= nu.sin(b)
    R[2,2]= nu.cos(b)
    invxyz= nu.array([[vx,vx,vx],
                    [vy,vy,vy],
                    [vz,vz,vz]])
    vrvlvb= (R.T*invxyz.T).sum(-1)
    vrvlvb[:,1]/= d*_K
    vrvlvb[:,2]/= d*_K
    return vrvlvb

@scalarDecorator
@degreeDecorator([],[0,1])
def XYZ_to_lbd(X,Y,Z,degree=False):
    """
    NAME:

       XYZ_to_lbd

    PURPOSE:

       transform from rectangular Galactic coordinates to spherical Galactic coordinates (works with vector inputs)

    INPUT:

       X - component towards the Galactic Center (in kpc; though this obviously does not matter))

       Y - component in the direction of Galactic rotation (in kpc)

       Z - component towards the North Galactic Pole (kpc)

       degree - (Bool) if True, return l and b in degrees

    OUTPUT:

       [l,b,d] in (rad or degree,rad or degree,kpc)

       For vector inputs [:,3]

    HISTORY:

       2009-10-24 - Written - Bovy (NYU)

       2014-06-14 - Re-written w/ numpy functions for speed and w/ decorators for beauty - Bovy (IAS)

    """
    #Whether to use degrees and scalar input is handled by decorators
    d= nu.sqrt(X**2.+Y**2.+Z**2.)
    b=nu.arcsin(Z/d)
    cosl= X/d/nu.cos(b)
    sinl= Y/d/nu.cos(b)
    l= nu.arcsin(sinl)
    l[cosl < 0.]= nu.pi-l[cosl < 0.]
    l[(cosl >= 0.)*(sinl < 0.)]+= 2.*nu.pi
    out= nu.empty((len(l),3))
    out[:,0]= l
    out[:,1]= b
    out[:,2]= d
    return out

@scalarDecorator
@degreeDecorator([2,3],[])
def pmrapmdec_to_pmllpmbb(pmra,pmdec,ra,dec,degree=False,epoch=2000.0):
    """
    NAME:

       pmrapmdec_to_pmllpmbb

    PURPOSE:

       rotate proper motions in (ra,dec) into proper motions in (l,b)

    INPUT:

       pmra - proper motion in ra (multplied with cos(dec)) [mas/yr]

       pmdec - proper motion in dec [mas/yr]

       ra - right ascension

       dec - declination

       degree - if True, ra and dec are given in degrees (default=False)

       epoch - epoch of ra,dec (right now only 2000.0 and 1950.0 are supported)

    OUTPUT:

       (pmll,pmbb) for vector inputs [:,2]

    HISTORY:

       2010-04-07 - Written - Bovy (NYU)

       2014-06-14 - Re-written w/ numpy functions for speed and w/ decorators for beauty - Bovy (IAS)

    """
    theta,dec_ngp,ra_ngp= get_epoch_angles(epoch)
    #Whether to use degrees and scalar input is handled by decorators
    dec[dec == dec_ngp]+= 10.**-16 #deal w/ pole.
    sindec_ngp= nu.sin(dec_ngp)
    cosdec_ngp= nu.cos(dec_ngp)
    sindec= nu.sin(dec)
    cosdec= nu.cos(dec)
    sinrarangp= nu.sin(ra-ra_ngp)
    cosrarangp= nu.cos(ra-ra_ngp)
    #These were replaced by Poleski (2013)'s equivalent form that is better at the poles
    #cosphi= (sindec_ngp-sindec*sinb)/cosdec/cosb
    #sinphi= sinrarangp*cosdec_ngp/cosb
    cosphi= sindec_ngp*cosdec-cosdec_ngp*sindec*cosrarangp
    sinphi= sinrarangp*cosdec_ngp
    norm= nu.sqrt(cosphi**2.+sinphi**2.)
    cosphi/= norm
    sinphi/= norm
    return (nu.array([[cosphi,-sinphi],[sinphi,cosphi]]).T\
                *nu.array([[pmra,pmra],[pmdec,pmdec]]).T).sum(-1)

@scalarDecorator
@degreeDecorator([2,3],[])
def pmllpmbb_to_pmrapmdec(pmll,pmbb,l,b,degree=False,epoch=2000.0):
    """
    NAME:

       pmllpmbb_to_pmrapmdec

    PURPOSE:

       rotate proper motions in (l,b) into proper motions in (ra,dec)

    INPUT:

       pmll - proper motion in l (multplied with cos(b)) [mas/yr]

       pmbb - proper motion in b [mas/yr]

       l - Galactic longitude

       b - Galactic lattitude

       degree - if True, l and b are given in degrees (default=False)

       epoch - epoch of ra,dec (right now only 2000.0 and 1950.0 are supported)

    OUTPUT:

       (pmra,pmdec), for vector inputs [:,2]

    HISTORY:

       2010-04-07 - Written - Bovy (NYU)

       2014-06-14 - Re-written w/ numpy functions for speed and w/ decorators for beauty - Bovy (IAS)

    """
    theta,dec_ngp,ra_ngp= get_epoch_angles(epoch)
    #Whether to use degrees and scalar input is handled by decorators
    radec = lb_to_radec(l,b,degree=False,epoch=epoch)
    ra= radec[:,0]
    dec= radec[:,1]
    dec[dec == dec_ngp]+= 10.**-16 #deal w/ pole.
    sindec_ngp= nu.sin(dec_ngp)
    cosdec_ngp= nu.cos(dec_ngp)
    sindec= nu.sin(dec)
    cosdec= nu.cos(dec)
    sinrarangp= nu.sin(ra-ra_ngp)
    cosrarangp= nu.cos(ra-ra_ngp)
    #These were replaced by Poleski (2013)'s equivalent form that is better at the poles
    #cosphi= (sindec_ngp-sindec*sinb)/cosdec/cosb
    #sinphi= sinrarangp*cosdec_ngp/cosb
    cosphi= sindec_ngp*cosdec-cosdec_ngp*sindec*cosrarangp
    sinphi= sinrarangp*cosdec_ngp
    norm= nu.sqrt(cosphi**2.+sinphi**2.)
    cosphi/= norm
    sinphi/= norm
    return (nu.array([[cosphi,sinphi],[-sinphi,cosphi]]).T\
                *nu.array([[pmll,pmll],[pmbb,pmbb]]).T).sum(-1)

def cov_pmrapmdec_to_pmllpmbb(cov_pmradec,ra,dec,degree=False,epoch=2000.0):
    """
    NAME:

       cov_pmrapmdec_to_pmllpmbb

    PURPOSE:

       propagate the proper motions errors through the rotation from (ra,dec) to (l,b)

    INPUT:

       covar_pmradec - uncertainty covariance matrix of the proper motion in ra (multplied with cos(dec)) and dec [2,2] or [:,2,2]

       ra - right ascension

       dec - declination

       degree - if True, ra and dec are given in degrees (default=False)

       epoch - epoch of ra,dec (right now only 2000.0 and 1950.0 are supported)

    OUTPUT:

       covar_pmllbb [2,2] or [:,2,2]

    HISTORY:

       2010-04-12 - Written - Bovy (NYU)

    """
    if len(cov_pmradec.shape) == 3:
        out= sc.zeros(cov_pmradec.shape)
        ndata= out.shape[0]
        lb = radec_to_lb(ra,dec,degree=degree,epoch=epoch)
        for ii in range(ndata):
            out[ii,:,:]= cov_pmradec_to_pmllbb_single(cov_pmradec[ii,:,:],
                                                      ra[ii],dec[ii],lb[ii,1],
                                                      degree,epoch)
        return out
    else:
        l,b = radec_to_lb(ra,dec,degree=degree,epoch=epoch)
        return cov_pmradec_to_pmllbb_single(cov_pmradec,ra,dec,b,degree,epoch)

def cov_pmradec_to_pmllbb_single(cov_pmradec,ra,dec,b,degree=False,epoch=2000.0):
    """
    NAME:
       cov_pmradec_to_pmllbb_single
    PURPOSE:
       propagate the proper motions errors through the rotation from (ra,dec)
       to (l,b) for scalar inputs
    INPUT:
       covar_pmradec - uncertainty covariance matrix of the proper motion
                      in ra (multplied with cos(dec)) and dec [2,2] or [:,2,2]
       ra - right ascension
       dec - declination
       degree - if True, ra and dec are given in degrees (default=False)
       epoch - epoch of ra,dec (right now only 2000.0 and 1950.0 are supported)
    OUTPUT:
       cov_pmllbb
    HISTORY:
       2010-04-12 - Written - Bovy (NYU)
    """
    theta,dec_ngp,ra_ngp= get_epoch_angles(epoch)
    if degree:
        sindec_ngp= m.sin(dec_ngp)
        cosdec_ngp= m.cos(dec_ngp)
        sindec= m.sin(dec*_DEGTORAD)
        cosdec= m.cos(dec*_DEGTORAD)
        sinrarangp= m.sin(ra*_DEGTORAD-ra_ngp)
        cosrarangp= m.cos(ra*_DEGTORAD-ra_ngp)
    else:
        sindec_ngp= m.sin(dec_ngp)
        cosdec_ngp= m.cos(dec_ngp)
        sindec= m.sin(dec)
        cosdec= m.cos(dec)
        sinrarangp= m.sin(ra-ra_ngp)
        cosrarangp= m.cos(ra-ra_ngp)
    #These were replaced by Poleski (2013)'s equivalent form that is better at the poles
    #cosphi= (sindec_ngp-sindec*sinb)/cosdec/cosb
    #sinphi= sinrarangp*cosdec_ngp/cosb
    cosphi= sindec_ngp*cosdec-cosdec_ngp*sindec*cosrarangp
    sinphi= sinrarangp*cosdec_ngp
    norm= m.sqrt(cosphi**2.+sinphi**2.)
    cosphi/= norm
    sinphi/= norm
    P= sc.array([[cosphi,sinphi],[-sinphi,cosphi]])
    return sc.dot(P,sc.dot(cov_pmradec,P.T))

def cov_dvrpmllbb_to_vxyz(d,e_d,e_vr,pmll,pmbb,cov_pmllbb,l,b,
                          plx=False,degree=False):
    """
    NAME:

       cov_dvrpmllbb_to_vxyz

    PURPOSE:

       propagate distance, radial velocity, and proper motion uncertainties to Galactic coordinates

    INPUT:

       d - distance [kpc, as/mas for plx]

       e_d - distance uncertainty [kpc, [as/mas] for plx]

       e_vr  - low velocity uncertainty [km/s]

       pmll - proper motion in l (*cos(b)) [ [as/mas]/yr ]

       pmbb - proper motion in b [ [as/mas]/yr ]

       cov_pmllbb - uncertainty covariance for proper motion

       l - Galactic longitude

       b - Galactic lattitude

    KEYWORDS:

       plx - if True, d is a parallax, and e_d is a parallax uncertainty

       degree - if True, l and b are given in degree

    OUTPUT:

       cov(vx,vy,vz) [3,3] or [:,3,3]

    HISTORY:

       2010-04-12 - Written - Bovy (NYU)

    """
    if plx:
        d= 1./d
        e_d*= d**2.
    if degree:
        l*= _DEGTORAD
        b*= _DEGTORAD
    if sc.array(d).shape == ():
        return cov_dvrpmllbb_to_vxyz_single(d,e_d,e_vr,pmll,pmbb,cov_pmllbb,
                                            l,b)
    else:
        ndata= len(d)
        out= sc.zeros((ndata,3,3))
        for ii in range(ndata):
            out[ii,:,:]= cov_dvrpmllbb_to_vxyz_single(d[ii],e_d[ii],e_vr[ii],
                                                      pmll[ii],pmbb[ii],
                                                      cov_pmllbb[ii,:,:],
                                                      l[ii],b[ii])

        return out
    
def cov_dvrpmllbb_to_vxyz_single(d,e_d,e_vr,pmll,pmbb,cov_pmllbb,l,b):
    """
    NAME:
       cov_dvrpmllbb_to_vxyz
    PURPOSE:
       propagate distance, radial velocity, and proper motion uncertainties to
       Galactic coordinates for scalar inputs
    INPUT:
       d - distance [kpc, as/mas for plx]
       e_d - distance uncertainty [kpc, [as/mas] for plx]
       e_vr  - low velocity uncertainty [km/s]
       pmll - proper motion in l (*cos(b)) [ [as/mas]/yr ]
       pmbb - proper motion in b [ [as/mas]/yr ]
       cov_pmllbb - uncertainty covariance for proper motion
       l - Galactic longitude [rad]
       b - Galactic lattitude [rad]
    OUTPUT:
       cov(vx,vy,vz) [3,3]
    HISTORY:
       2010-04-12 - Written - Bovy (NYU)
    """
    M= _K*sc.array([[pmll,d,0.],[pmbb,0.,d]])
    cov_dpmllbb= sc.zeros((3,3))
    cov_dpmllbb[0,0]= e_d**2.
    cov_dpmllbb[1:3,1:3]= cov_pmllbb
    cov_vlvb= sc.dot(M,sc.dot(cov_dpmllbb,M.T))
    cov_vrvlvb= sc.zeros((3,3))
    cov_vrvlvb[0,0]= e_vr**2.
    cov_vrvlvb[1:3,1:3]= cov_vlvb
    R= sc.array([[m.cos(l)*m.cos(b), m.sin(l)*m.cos(b), m.sin(b)],
                 [-m.sin(l),m.cos(l),0.],
                 [-m.cos(l)*m.sin(b),-m.sin(l)*m.sin(b), m.cos(b)]])
    return sc.dot(R.T,sc.dot(cov_vrvlvb,R))

def XYZ_to_galcenrect(X,Y,Z,Xsun=1.,Ysun=0.,Zsun=0.):
    """
    NAME:

       XYZ_to_galcenrect

    PURPOSE:

       transform XYZ coordinates (wrt Sun) to rectangular Galactocentric coordinates

    INPUT:

       X - X

       Y - Y

       Z - Z

    OUTPUT:

       (Xg, Yg, Zg)

    HISTORY:

       2010-09-24 - Written - Bovy (NYU)

    """
    return (-X+Xsun,Y+Ysun,Z+Zsun)

def galcenrect_to_XYZ(X,Y,Z,Xsun=1.,Ysun=0.,Zsun=0.):
    """
    NAME:

       galcenrect_to_XYZ

    PURPOSE:

       transform rectangular Galactocentric to XYZ coordinates (wrt Sun) coordinates

    INPUT:

       X, Y, Z - Galactocentric rectangular coordinates

    OUTPUT:

       (X, Y, Z)

    HISTORY:

       2011-02-23 - Written - Bovy (NYU)

    """
    return (-X+Xsun,Y-Ysun,Z-Zsun)

def rect_to_cyl(X,Y,Z):
    """
    NAME:

       rect_to_cyl

    PURPOSE:

       convert from rectangular to cylindrical coordinates

    INPUT:

       X, Y, Z - rectangular coordinates

    OUTPUT:

       [:,3] R,phi,z

    HISTORY:

       2010-09-24 - Written - Bovy (NYU)

    """
    R= sc.sqrt(X**2.+Y**2.)
    phi= sc.arcsin(Y/R)
    if isinstance(X,float) and X < 0.:
        phi= m.pi-phi
    elif isinstance(X,sc.ndarray):
        phi[(X < 0.)]= m.pi-phi[(X < 0.)]
    return (R,phi,Z)

def cyl_to_rect(R,phi,Z):
    """
    NAME:

       cyl_to_rect

    PURPOSE:

       convert from cylindrical to rectangular coordinates

    INPUT:

       R, phi, Z - cylindrical coordinates

    OUTPUT:

       [:,3] X,Y,Z

    HISTORY:

       2011-02-23 - Written - Bovy (NYU)

    """
    return (R*sc.cos(phi),R*sc.sin(phi),Z)

def XYZ_to_galcencyl(X,Y,Z,Xsun=1.,Ysun=0.,Zsun=0.):
    """
    NAME:

       XYZ_to_galcencyl

    PURPOSE:

       transform XYZ coordinates (wrt Sun) to cylindrical Galactocentric coordinates

    INPUT:

       X - X

       Y - Y

       Z - Z

    OUTPUT:

       [:,3]= R,phi,z

    HISTORY:

       2010-09-24 - Written - Bovy (NYU)

    """
    Xg,Yg,Zg= XYZ_to_galcenrect(X,Y,Z,Xsun=Xsun,Ysun=Ysun,Zsun=Zsun)
    return rect_to_cyl(Xg,Yg,Zg)
    
def galcencyl_to_XYZ(R,phi,Z,Xsun=1.,Ysun=0.,Zsun=0.):
    """
    NAME:

       galcencyl_to_XYZ

    PURPOSE:

       transform cylindrical Galactocentric coordinates to XYZ coordinates (wrt Sun)

    INPUT:

       R, phi, Z - Galactocentric cylindrical coordinates

    OUTPUT:

       [:,3]= X,Y,Z

    HISTORY:

       2011-02-23 - Written - Bovy (NYU)

    """
    Xr,Yr,Zr= cyl_to_rect(R,phi,Z)
    return galcenrect_to_XYZ(Xr,Yr,Zr,Xsun=Xsun,Ysun=Ysun,Zsun=Zsun)
    
def vxvyvz_to_galcenrect(vx,vy,vz,vsun=[0.,1.,0.]):
    """
    NAME:

       vxvyvz_to_galcenrect

    PURPOSE:

       transform velocities in XYZ coordinates (wrt Sun) to rectangular Galactocentric coordinates for velocities

    INPUT:

       vx - U

       vy - V

       vz - W

       vsun - velocity of the sun in the GC frame ndarray[3]

    OUTPUT:

       [:,3]= vXg, vYg, vZg

    HISTORY:

       2010-09-24 - Written - Bovy (NYU)

    """
    return sc.array([-vx+vsun[0],vy+vsun[1],vz+vsun[2]])

def vxvyvz_to_galcencyl(vx,vy,vz,X,Y,Z,vsun=[0.,1.,0.],galcen=False):
    """
    NAME:

       vxvyvz_to_galcencyl

    PURPOSE:

       transform velocities in XYZ coordinates (wrt Sun) to cylindrical Galactocentric coordinates for velocities

    INPUT:

       vx - U

       vy - V

       vz - W

       X - X in Galactocentric rectangular coordinates

       Y - Y in Galactocentric rectangular coordinates

       Z - Z in Galactocentric rectangular coordinates

       vsun - velocity of the sun in the GC frame ndarray[3]

       galcen - if True, then X,Y,Z are in cylindrical Galactocentric coordinates rather than rectangular coordinates

    OUTPUT:

       vRg, vTg, vZg

    HISTORY:

       2010-09-24 - Written - Bovy (NYU)

    """
    vx,vy,vz= vxvyvz_to_galcenrect(vx,vy,vz,vsun=vsun)
    return rect_to_cyl_vec(vx,vy,vz,X,Y,Z,cyl=galcen)

def galcenrect_to_vxvyvz(vXg,vYg,vZg,vsun=[0.,1.,0.]):
    """
    NAME:

       galcenrect_to_vxvyvz

    PURPOSE:

       transform rectangular Galactocentric coordinates to XYZ coordinates (wrt Sun) for velocities

    INPUT:

       vXg - Galactocentric x-velocity

       vYg - Galactocentric y-velocity

       vZg - Galactocentric z-velocity

       vsun - velocity of the sun in the GC frame ndarray[3]

    OUTPUT:

       [:,3]= vx, vy, vz

    HISTORY:

       2011-02-24 - Written - Bovy (NYU)

    """
    return sc.array([-vXg+vsun[0],vYg-vsun[1],vZg-vsun[2]])
#Old form, does not seem to be necessary anymore? Bovy 06/08/14
"""
    try:
        return sc.array([-vXg+vsun[0],vYg-vsun[1],vZg-vsun[2]])
    except ValueError: #annoying bug for one-d, make sure they are arrays
        return sc.array([-sc.array([vXg]).flatten()[0]+vsun[0],
                          sc.array([vYg]).flatten()[0]-vsun[1],
                          sc.array([vZg]).flatten()[0]-vsun[2]])      
"""

def galcencyl_to_vxvyvz(vR,vT,vZ,phi,vsun=[0.,1.,0.]):
    """
    NAME:

       galcencyl_to_vxvyvz

    PURPOSE:

       transform cylindrical Galactocentric coordinates to XYZ (wrt Sun) coordinates for velocities

    INPUT:

       vR - Galactocentric radial velocity

       vT - Galactocentric tangential velocity

       vZ - Galactocentric vertical velocity

       phi - Galactocentric azimuth

       vsun - velocity of the sun in the GC frame ndarray[3]

    OUTPUT:

       vx,vy,vz

    HISTORY:

       2011-02-24 - Written - Bovy (NYU)

    """
    vXg, vYg, vZg= cyl_to_rect_vec(vR,vT,vZ,phi)
    return galcenrect_to_vxvyvz(vXg,vYg,vZg,vsun=vsun)

def rect_to_cyl_vec(vx,vy,vz,X,Y,Z,cyl=False):
    """
    NAME:

       rect_to_cyl_vec

    PURPOSE:

       transform vectors from rectangular to cylindrical coordinates vectors

    INPUT:

       vx - 

       vy - 

       vz - 

       X - X

       Y - Y

       Z - Z

       cyl - if True, X,Y,Z are already cylindrical

    OUTPUT:

       vR,vT,vz

    HISTORY:

       2010-09-24 - Written - Bovy (NYU)

    """
    if not cyl:
        R,phi,Z= rect_to_cyl(X,Y,Z)
    else:
        phi= Y
    vr=+vx*sc.cos(phi)+vy*sc.sin(phi)
    vt= -vx*sc.sin(phi)+vy*sc.cos(phi)
    return (vr,vt,vz)

def cyl_to_rect_vec(vr,vt,vz,phi):
    """
    NAME:

       cyl_to_rect_vec

    PURPOSE:

       transform vectors from cylindrical to rectangular coordinate vectors

    INPUT:

       vr - radial velocity

       vt - tangential velocity

       vz - vertical velocity

       phi - azimuth

    OUTPUT:

       vx,vy,vz

    HISTORY:

       2011-02-24 - Written - Bovy (NYU)

    """
    vx= vr*sc.cos(phi)-vt*sc.sin(phi)
    vy= vr*sc.sin(phi)+vt*sc.cos(phi)
    return (vx,vy,vz)

def cyl_to_rect_jac(*args):
    """
    NAME:

       cyl_to_rect_jac

    PURPOSE:

       calculate the Jacobian of the cylindrical to rectangular conversion

    INPUT:

       R, phi, Z- cylindrical coordinates

       vR, vT, vZ- cylindrical velocities

       if 6 inputs: R,vR,vT,z,vz,phi

       if 3: R, phi, Z

    OUTPUT:

       jacobian d(rect)/d(cyl)

    HISTORY:

       2013-12-09 - Written - Bovy (IAS)

    """
    out= sc.zeros((6,6))
    if len(args) == 3:
        R, phi, Z= args
        vR, vT, vZ= 0., 0., 0.
        outIndx= sc.array([True,False,False,True,False,True],dtype='bool')
    elif len(args) == 6:
        R, vR, vT, Z, vZ, phi= args
        outIndx= sc.ones(6,dtype='bool')
    cp= sc.cos(phi)
    sp= sc.sin(phi)
    out[0,0]= cp
    out[0,5]= -R*sp
    out[1,0]= sp
    out[1,5]= R*cp
    out[2,3]= 1.
    out[3,1]= cp
    out[3,2]= -sp
    out[3,5]= -vT*cp-vR*sp
    out[4,1]= sp
    out[4,2]= cp
    out[4,5]= -vT*sp+vR*cp
    out[5,4]= 1.
    if len(args) == 3:
        out= out[:3,outIndx]
        out[:,[1,2]]= out[:,[2,1]]
    return out

def galcenrect_to_XYZ_jac(*args,**kwargs):
    """
    NAME:

       galcenrect_to_XYZ_jac
    PURPOSE:

       calculate the Jacobian of the Galactocentric rectangular to Galactic coordinates

    INPUT:

       X,Y,Z- Galactocentric rectangular coordinates

       vX, vY, vZ- Galactocentric rectangular velocities

       if 6 inputs: X,Y,Z,vX,vY,vZ

       if 3: X,Y,Z

    OUTPUT:

       jacobian d(galcen.)/d(Galactic)

    HISTORY:

       2013-12-09 - Written - Bovy (IAS)

    """
    out= sc.zeros((6,6))
    out[0,0]= -1.
    out[1,1]= 1.
    out[2,2]= 1.
    if len(args) == 3: return out[:3,:3]
    out[3,3]= -1.
    out[4,4]= 1.
    out[5,5]= 1.
    return out

def lbd_to_XYZ_jac(*args,**kwargs):
    """
    NAME:

       lbd_to_XYZ_jac

    PURPOSE:

       calculate the Jacobian of the Galactic spherical coordinates to Galactic rectangular coordinates transformation

    INPUT:

       l,b,D- Galactic spherical coordinates

       vlos,pmll,pmbb- Galactic spherical velocities (some as proper motions)

       if 6 inputs: l,b,D,vlos,pmll,pmbb

       if 3: l,b,D

       degree= (False) if True, l and b are in degrees

    OUTPUT:

       jacobian

    HISTORY:

       2013-12-09 - Written - Bovy (IAS)

    """
    out= sc.zeros((6,6))
    if len(args) == 3:
        l,b,D= args
        vlos, pmll, pmbb= 0., 0., 0.
    elif len(args) == 6:
        l,b,D,vlos,pmll,pmbb= args
    if kwargs.has_key('degree') and kwargs['degree']:
        l*= _DEGTORAD
        b*= _DEGTORAD
    cl= sc.cos(l)
    sl= sc.sin(l)
    cb= sc.cos(b)
    sb= sc.sin(b)
    out[0,0]= -D*cb*sl
    out[0,1]= -D*sb*cl
    out[0,2]= cb*cl
    out[1,0]= D*cb*cl
    out[1,1]= -D*sb*sl
    out[1,2]= cb*sl
    out[2,1]= D*cb
    out[2,2]= sb
    if len(args) == 3:
        if kwargs.has_key('degree') and kwargs['degree']:
            out[:,0]*= _DEGTORAD
            out[:,1]*= _DEGTORAD
        return out[:3,:3]
    out[3,0]= -sl*cb*vlos-cl*_K*D*pmll+sb*sl*_K*D*pmbb
    out[3,1]= -cl*sb*vlos-cb*cl*_K*D*pmbb
    out[3,2]= -sl*_K*pmll-sb*cl*_K*pmbb
    out[3,3]= cl*cb
    out[3,4]= -sl*_K*D
    out[3,5]= -cl*sb*_K*D
    out[4,0]= cl*cb*vlos-sl*_K*D*pmll-cl*sb*_K*D*pmbb
    out[4,1]= -sl*sb*vlos-sl*cb*_K*D*pmbb
    out[4,2]= cl*_K*pmll-sl*sb*_K*pmbb
    out[4,3]= sl*cb
    out[4,4]= cl*_K*D
    out[4,5]= -sl*sb*_K*D
    out[5,1]= cb*vlos-sb*_K*D*pmbb
    out[5,2]= cb*_K*pmbb
    out[5,3]= sb
    out[5,5]= cb*_K*D
    if kwargs.has_key('degree') and kwargs['degree']:
        out[:,0]*= _DEGTORAD
        out[:,1]*= _DEGTORAD
    return out

def dl_to_rphi_2d(d,l,degree=False,ro=1.,phio=0.):
    """
    NAME:

       dl_to_rphi_2d

    PURPOSE:

       convert Galactic longitude and distance to Galactocentric radius and azimuth

    INPUT:

       d - distance

       l - Galactic longitude [rad/deg if degree]

    KEYWORDS:

       degree= (False): l is in degrees rather than rad

       ro= (1) Galactocentric radius of the observer

       phio= (0) Galactocentric azimuth of the observer [rad/deg if degree]

    OUTPUT:

       (R,phi); phi in degree if degree

    HISTORY:

       2012-01-04 - Written - Bovy (IAS)

    """
    scalarOut, listOut= False, False
    if isinstance(d,(int,float)):
        d= sc.array([d])
        scalarOut= True
    elif isinstance(d,list):
        d= sc.array(d)
        listOut= True
    if isinstance(l,(int,float)):
        l= sc.array([l])
    elif isinstance(l,list):
        l= sc.array(l)
    if degree:
        l*= _DEGTORAD
    R= sc.sqrt(ro**2.+d**2.-2.*d*ro*sc.cos(l))
    phi= sc.arcsin(d/R*sc.sin(l))
    indx= (ro/sc.cos(l) < d)*(sc.cos(l) > 0.)
    phi[indx]= sc.pi-sc.arcsin(d[indx]/R[indx]*sc.sin(l[indx]))
    if degree:
        phi/= _DEGTORAD
    phi+= phio
    if scalarOut:
        return (R[0],phi[0])
    elif listOut:
        return (list(R),list(phi))
    else:
        return (R,phi)

def rphi_to_dl_2d(R,phi,degree=False,ro=1.,phio=0.):
    """
    NAME:

       rphi_to_dl_2d

    PURPOSE:

       convert Galactocentric radius and azimuth to distance and Galactic longitude

    INPUT:

       R - Galactocentric radius

       phi - Galactocentric azimuth [rad/deg if degree]

    KEYWORDS:

       degree= (False): phi is in degrees rather than rad

       ro= (1) Galactocentric radius of the observer

       phio= (0) Galactocentric azimuth of the observer [rad/deg if degree]

    OUTPUT:

       (d,l); phi in degree if degree

    HISTORY:

       2012-01-04 - Written - Bovy (IAS)

    """
    scalarOut, listOut= False, False
    if isinstance(R,(int,float)):
        R= sc.array([R])
        scalarOut= True
    elif isinstance(R,list):
        R= sc.array(R)
        listOut= True
    if isinstance(phi,(int,float)):
        phi= sc.array([phi])
    elif isinstance(phi,list):
        phi= sc.array(phi)
    phi-= phio
    if degree:
        phi*= _DEGTORAD
    d= sc.sqrt(R**2.+ro**2.-2.*R*ro*sc.cos(phi))
    l= sc.arcsin(R/d*sc.sin(phi))
    indx= (ro/sc.cos(phi) < R)*(sc.cos(phi) > 0.)
    l[indx]= sc.pi-sc.arcsin(R[indx]/d[indx]*sc.sin(phi[indx]))
    if degree:
        l/= _DEGTORAD
    if scalarOut:
        return (d[0],l[0])
    elif listOut:
        return (list(d),list(l))
    else:
        return (d,l)

def Rz_to_coshucosv(R,z,delta=1.):
    """
    NAME:

       Rz_to_coshucosv

    PURPOSE:

       calculate prolate confocal cosh(u) and cos(v) coordinates from R,z, and delta

    INPUT:

       R - radius

       z - height

       delta= focus

    OUTPUT:

       (cosh(u),cos(v))

    HISTORY:

       2012-11-27 - Written - Bovy (IAS)

    """
    d12= (z+delta)**2.+R**2.
    d22= (z-delta)**2.+R**2.
    coshu= 0.5/delta*(sc.sqrt(d12)+sc.sqrt(d22))
    cosv=  0.5/delta*(sc.sqrt(d12)-sc.sqrt(d22))
    return (coshu,cosv)

def Rz_to_uv(R,z,delta=1.):
    """
    NAME:

       Rz_to_uv

    PURPOSE:

       calculate prolate confocal u and v coordinates from R,z, and delta

    INPUT:

       R - radius

       z - height

       delta= focus

    OUTPUT:

       (u,v)

    HISTORY:

       2012-11-27 - Written - Bovy (IAS)

    """
    coshu, cosv= Rz_to_coshucosv(R,z,delta)
    u= sc.arccosh(coshu)
    v= sc.arccos(cosv)
    return (u,v)

def uv_to_Rz(u,v,delta=1.):
    """
    NAME:

       uv_to_Rz

    PURPOSE:

       calculate R and z from prolate confocal u and v coordinates

    INPUT:

       u - confocal u

       v - confocal v

       delta= focus

    OUTPUT:

       (R,z)

    HISTORY:

       2012-11-27 - Written - Bovy (IAS)

    """
    R= delta*sc.sinh(u)*sc.sin(v)
    z= delta*sc.cosh(u)*sc.cos(v)
    return (R,z)

def get_epoch_angles(epoch=2000.0):
    """
    NAME:

       get_epoch_angles

    PURPOSE:

       get the angles relevant for the transformation from ra, dec to l,b for the given epoch
    INPUT:

       epoch - epoch of ra,dec (right now only 2000.0 and 1950.0 are supported

    OUTPUT:

       set of angles

    HISTORY:

       2010-04-07 - Written - Bovy (NYU)

    """
    if epoch == 2000.0:
        theta= 122.932/180.*sc.pi
        dec_ngp= 27.12825/180.*sc.pi
        ra_ngp= 192.85948/180.*sc.pi
    elif epoch == 1950.0:
        theta= 123./180.*sc.pi
        dec_ngp= 27.4/180.*sc.pi
        ra_ngp= 192.25/180.*sc.pi
    else:
        raise IOError("Only epochs 1950 and 2000 are supported")
    return (theta,dec_ngp,ra_ngp)

