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
#
##############################################################################
#############################################################################
#Copyright (c) 2010 - 2011, Jo Bovy
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
import math as m
import scipy as sc
_DEGTORAD= m.pi/180.
_K=4.74047
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
    """
    #First calculate the transformation matrix T
    theta,dec_ngp,ra_ngp= get_epoch_angles(epoch)
    T= sc.dot(sc.array([[sc.cos(theta),sc.sin(theta),0.],[sc.sin(theta),-sc.cos(theta),0.],[0.,0.,1.]]),sc.dot(sc.array([[-sc.sin(dec_ngp),0.,sc.cos(dec_ngp)],[0.,1.,0.],[sc.cos(dec_ngp),0.,sc.sin(dec_ngp)]]),sc.array([[sc.cos(ra_ngp),sc.sin(ra_ngp),0.],[-sc.sin(ra_ngp),sc.cos(ra_ngp),0.],[0.,0.,1.]])))
    transform={}
    transform['T']= T
    if sc.array(ra).shape == ():
        return radec_to_lb_single(ra,dec,transform,degree)
    else:
        function= sc.frompyfunc(radec_to_lb_single,4,2)
        return sc.array(function(ra,dec,transform,degree)).T

def radec_to_lb_single(ra,dec,T,degree=False):
    """
    NAME:
       radec_to_lb_single
    PURPOSE:
       transform from equatorial coordinates to Galactic coordinates for a single pair of ra,dec
    INPUT:
       ra - right ascension
       dec - declination
       T - epoch dependent transformation matrix (dictionary)
       degree - (Bool) if True, ra and dec are given in degree and l and b will be as well
    OUTPUT:
       l,b
    HISTORY:
       2009-11-12 - Written - Bovy (NYU)
    """
    T=T['T']
    if degree:
        thisra= ra/180.*sc.pi
        thisdec= dec/180.*sc.pi
    else:
        thisra= ra
        thisdec= dec
    XYZ=sc.array([sc.cos(thisdec)*sc.cos(thisra),sc.cos(thisdec)*sc.sin(thisra),sc.sin(thisdec)])
    galXYZ= sc.dot(T,XYZ)
    b= m.asin(galXYZ[2])
    l= m.atan(galXYZ[1]/galXYZ[0])
    if galXYZ[0]/sc.cos(b) < 0.:
        l+= sc.pi
    if l < 0.:
        l+= 2.*sc.pi
    if degree:
        return (l/sc.pi*180.,b/sc.pi*180.)
    else:
        return (l,b)

def lb_to_radec(l,b,degree=False,epoch=2000.0):
    """
    NAME:
       lb_to_radec
    PURPOSE:
       transform from Galactic coordinates to equatorial coordinates
    INPUT:
       l - Galactic longitude
       b - Galactic lattitude
       degree - (Bool) if True, l and b are given in degree and ra and dec
                will be as well
       epoch - epoch of target ra,dec
               (right now only 2000.0 and 1950.0 are supported)
    OUTPUT:
       ra,dec
       For vector inputs [:,2]
    HISTORY:
       2010-04-07 - Written - Bovy (NYU)
    """
    #First calculate the transformation matrix T'
    theta,dec_ngp,ra_ngp= get_epoch_angles(epoch)
    T= sc.dot(sc.array([[sc.cos(ra_ngp),-sc.sin(ra_ngp),0.],[sc.sin(ra_ngp),sc.cos(ra_ngp),0.],[0.,0.,1.]]),sc.dot(sc.array([[-sc.sin(dec_ngp),0.,sc.cos(dec_ngp)],[0.,1.,0.],[sc.cos(dec_ngp),0.,sc.sin(dec_ngp)]]),sc.array([[sc.cos(theta),sc.sin(theta),0.],[sc.sin(theta),-sc.cos(theta),0.],[0.,0.,1.]])))
    transform={}
    transform['T']= T
    if sc.array(l).shape == ():
        return lb_to_radec_single(l,b,transform,degree)
    else:
        function= sc.frompyfunc(lb_to_radec_single,4,2)
        return sc.array(function(l,b,transform,degree)).T

def lb_to_radec_single(l,b,T,degree=False):
    """
    NAME:
       lb_to_radec_single
    PURPOSE:
       transform from Galactic coordinates to equatorial coordinates
       for a single pair of l,b
    INPUT:
       l - Galactic longitude
       b - Galactic lattitude
       T - epoch dependent transformation matrix (dictionary)
       degree - (Bool) if True, l and b are given in degree and ra and dec
                 will be as well
    OUTPUT:
       ra,dec
    HISTORY:
       2010-04-08 - Written - Bovy (NYU)
    """
    T=T['T']
    if degree:
        thisl= l*_DEGTORAD
        thisb= b*_DEGTORAD
    else:
        thisl= l
        thisb= b
    XYZ=sc.array([m.cos(thisb)*m.cos(thisl),m.cos(thisb)*m.sin(thisl),m.sin(thisb)])
    eqXYZ= sc.dot(T,XYZ)
    dec= m.asin(eqXYZ[2])
    ra= m.atan(eqXYZ[1]/eqXYZ[0])
    if eqXYZ[0] < 0.:
        ra+= m.pi
    elif ra < 0.:
        ra+= 2.*m.pi
    if degree:
        return (ra/m.pi*180.,dec/m.pi*180.)
    else:
        return (ra,dec)

def lbd_to_XYZ(l,b,d,degree=False):
    """
    NAME:
       lbd_to_XYZ
    PURPOSE:
       transform from spherical Galactic coordinates to rectangular Galactic coordinates
       works with vector inputs
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
    """
    if sc.array(l).shape == ():
        return lbd_to_XYZ_single(l,b,d,degree)
    else:
        function= sc.frompyfunc(lbd_to_XYZ_single,4,3)
        return sc.array(function(l,b,d,degree),dtype=sc.float64).T

def lbd_to_XYZ_single(l,b,d,degree=False):
    """
    NAME:
       lbd_to_XYZ_single
    PURPOSE:
       transform from spherical Galactic coordinates to rectangular Galactic coordinates
       works with vector inputs
    INPUT:
       l - Galactic longitude (rad)
       b - Galactic lattitude (rad)
       d - distance (arbitrary units)
       degree - (bool) if True, l and b are in degrees
    OUTPUT:
       [X,Y,Z] in whatever units d was in
    HISTORY:
       2009-10-24- Written - Bovy (NYU)
    """
    if degree:
        l= l*m.pi/180.
        b= b*m.pi/180.
    return (d*m.cos(b)*m.cos(l),d*m.cos(b)*m.sin(l),d*m.sin(b))

def rectgal_to_sphergal(X,Y,Z,vx,vy,vz,degree=False):
    """
    NAME:
       rectgal_to_sphergal
    PURPOSE:
       transform phase-space coordinates in rectangular Galactic coordinates to
       spherical Galactic coordinates
       Can take vector inputs
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
       transform phase-space coordinates in spherical Galactic coordinates to
       rectangular Galactic coordinates
       Can take vector inputs
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

def vrpmllpmbb_to_vxvyvz(vr,pmll,pmbb,l,b,d,XYZ=False,degree=False):
    """
    NAME:
       vrpmllpmbb_to_vxvyvz
    PURPOSE:
       Transform velocities in the spherical Galactic coordinate frame to the rectangular Galactic coordinate frame
       Can take vector inputs
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
    HISTORY:
       2009-10-24 - Written - Bovy (NYU)
    """
    if sc.array(l).shape == ():
        return vrpmllpmbb_to_vxvyvz_single(vr,pmll,pmbb,l,b,d,XYZ,degree)
    else:
        function= sc.frompyfunc(vrpmllpmbb_to_vxvyvz_single,8,3)
        return sc.array(function(vr,pmll,pmbb,l,b,d,XYZ,degree)).T
    
def vrpmllpmbb_to_vxvyvz_single(vr,pmll,pmbb,l,b,d,XYZ,degree):
    """
    NAME:
       vrpmllpmbb_to_vxvyvz
    PURPOSE:
       Transform velocities in the spherical Galactic coordinate frame to the rectangular Galactic coordinate frame
       Can take vector inputs
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
    HISTORY:
       2009-10-24 - Written - Bovy (NYU)
    """
    if XYZ:
        lbd= XYZ_to_lbd(l,b,d,degree=degree)
        if degree:
            l= lbd[0]*_DEGTORAD
            b= lbd[1]*_DEGTORAD
        else:
            l= lbd[0]
            b= lbd[1]
        d= lbd[2]
    else:
        if degree:
            l*= _DEGTORAD
            b*= _DEGTORAD
    R=sc.zeros((3,3))
    R[0,0]= m.cos(l)*m.cos(b)
    R[1,0]= -m.sin(l)
    R[2,0]= -m.cos(l)*m.sin(b)
    R[0,1]= m.sin(l)*m.cos(b)
    R[1,1]= m.cos(l)
    R[2,1]= -m.sin(l)*m.sin(b)
    R[0,2]= m.sin(b)
    R[2,2]= m.cos(b)
    vxvyvz= sc.dot(R.T,sc.array([vr,d*pmll*_K,d*pmbb*_K]))
    return (vxvyvz[0],vxvyvz[1],vxvyvz[2])

def vxvyvz_to_vrpmllpmbb(vx,vy,vz,l,b,d,XYZ=False,degree=False):
    """
    NAME:
       vxvyvz_to_vrpmllpmbb
    PURPOSE:
       Transform velocities in the rectangular Galactic coordinate frame to the spherical Galactic coordinate frame
       Can take vector inputs
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
    HISTORY:
       2009-10-24 - Written - Bovy (NYU)
    """
    if sc.array(l).shape == ():
        return vxvyvz_to_vrpmllpmbb_single(vx,vy,vz,l,b,d,XYZ,degree)
    else:
        function= sc.frompyfunc(vxvyvz_to_vrpmllpmbb_single,8,3)
        return sc.array(function(vx,vy,vz,l,b,d,XYZ,degree)).T
    
def vxvyvz_to_vrpmllpmbb_single(vx,vy,vz,l,b,d,XYZ=False,degree=False):
    """
    NAME:
       vxvyvz_to_vrpmllpmbb_single
    PURPOSE:
       Transform velocities in the rectangular Galactic coordinate frame to the spherical Galactic coordinate frame
    INPUT:
       vx - velocity towards the Galactic Center (km/s)
       vy - velocity in the direction of Galactic rotation (km/s)
       vz - velocity towards the North Galactic Pole (km/s)
       l - Galactic longitude
       b - Galactic lattitude
       d - distance (kpc)
       XYZ - (bool) If True, then l,b,d is actually X,Y,Z
             (rectangular Galactic coordinates)
       degree - (bool) if True, l and b are in degrees
    OUTPUT:
       (vr,pmll,pmbb) in (km/s,mas/yr,mas/yr); pmll = mu_l * cos(b)
    HISTORY:
       2009-10-24 - Written - Bovy (NYU)
    """ 
    if XYZ:
        lbd= XYZ_to_lbd(l,b,d,degree)
        if degree:
            l= lbd[0]*_DEGTORAD
            b= lbd[1]*_DEGTORAD
        else:
            l= lbd[0]
            b= lbd[1]
        d= lbd[2]
    else:
        if degree:
            l*= _DEGTORAD
            b*= _DEGTORAD
    R=sc.zeros((3,3))
    R[0,0]= m.cos(l)*m.cos(b)
    R[1,0]= -m.sin(l)
    R[2,0]= -m.cos(l)*m.sin(b)
    R[0,1]= m.sin(l)*m.cos(b)
    R[1,1]= m.cos(l)
    R[2,1]= -m.sin(l)*m.sin(b)
    R[0,2]= m.sin(b)
    R[2,2]= m.cos(b)
    vrvlvb= sc.dot(R,sc.array([vx,vy,vz]))
    pmll= vrvlvb[1]/d/_K
    pmbb= vrvlvb[2]/d/_K
    return (vrvlvb[0],pmll,pmbb)

def XYZ_to_lbd(X,Y,Z,degree=False):
    """
    NAME:
       XYZ_to_lbd
    PURPOSE:
       transform from rectangular Galactic coordinates to spherical Galactic coordinates
       works with vector inputs
    INPUT:
       X - component towards the Galactic Center (in kpc; though this obviously does not matter))
       Y - component in the direction of Galactic rotation (in kpc)
       Z - component towards the North Galactic Pole (kpc)
       degree - (Bool) if True, return l and b in degrees
    OUTPUT:
       [l,b,d] in (rad,rad,kpc)
       For vector inputs [:,3]
    HISTORY:
       2009-10-24 - Written - Bovy (NYU)
    """
    if sc.array(X).shape == ():
        return XYZ_to_lbd_single(X,Y,Z,degree)
    else:
        function= sc.frompyfunc(XYZ_to_lbd_single,4,3)
        return sc.array(function(X,Y,Z,degree)).T

def XYZ_to_lbd_single(X,Y,Z,degree):
    """
    NAME:
       XYZ_to_lbd_single
    PURPOSE:
       transform a single coordinate in rectangular Galactic coordinates to spherical Galactic coordinates
    INPUT:
       X - component towards the Galactic Center
       Y - component in the direction of Galactic rotation
       Z - component towards the North Galactic Pole
       degree - (Bool) if True, return l and b in degrees
    OUTPUT:
       [l,b,d] in (rad,rad,kpc)
    HISTORY:
       2009-10-24 - Written - Bovy (NYU)
    """
    d= sc.sqrt(X**2.+Y**2.+Z**2.)
    b=m.asin(Z/d)
    cosl= X/d/m.cos(b)
    sinl= Y/d/m.cos(b)
    l= m.asin(sinl)
    if cosl < 0:
        l= m.pi-l
    elif sinl < 0.:
        l=2.*m.pi+l
    if degree:
        return (l/m.pi*180.,b/m.pi*180.,d)
    else:
        return (l,b,d)

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
    """
    if sc.array(pmra).shape == ():
        l,b = radec_to_lb(ra,dec,degree=degree,epoch=epoch)
        return pmrapmdec_to_pmllpmbb_single(pmra,pmdec,ra,dec,b,degree,epoch)
    else:
        lb = radec_to_lb(ra,dec,degree=degree,epoch=epoch)
        function= sc.frompyfunc(pmrapmdec_to_pmllpmbb_single,7,2)
        return sc.array(function(pmra,pmdec,ra,dec,lb[:,1],degree,epoch)).T

def pmrapmdec_to_pmllpmbb_single(pmra,pmdec,ra,dec,b,degree=False,epoch=2000.0):
    """
    NAME:
       pmrapmdec_to_pmllpmbb_single
    PURPOSE:
       rotate proper motions in (ra,dec) into proper motions in (l,b) for
       scalar inputs
    INPUT:
       pmra - proper motion in ra (multplied with cos(dec)) [mas/yr]
       pmdec - proper motion in dec [mas/yr]
       ra - right ascension
       dec - declination
       b - Galactic lattitude
       degree - if True, ra and dec are given in degrees (default=False)
       epoch - epoch of ra,dec (right now only 2000.0 and 1950.0 are supported)
    OUTPUT:
       (pmll,pmbb)
    HISTORY:
       2010-04-07 - Written - Bovy (NYU)
    """
    theta,dec_ngp,ra_ngp= get_epoch_angles(epoch)
    if degree:
        sindec_ngp= m.sin(dec_ngp)
        cosdec_ngp= m.cos(dec_ngp)
        sindec= m.sin(dec*_DEGTORAD)
        sinb= m.sin(b*_DEGTORAD)
        cosdec= m.cos(dec*_DEGTORAD)
        cosb= m.cos(b*_DEGTORAD)
        sinrarangp= m.sin(ra*_DEGTORAD-ra_ngp)
    else:
        sindec_ngp= m.sin(dec_ngp)
        cosdec_ngp= m.cos(dec_ngp)
        sindec= m.sin(dec)
        sinb= m.sin(b)
        cosdec= m.cos(dec)
        cosb= m.sin(b)
        sinrarangp= m.sin(ra*-ra_ngp)
    cosphi= (sindec_ngp-sindec*sinb)/cosdec/cosb
    sinphi= sinrarangp*cosdec_ngp/cosb
    out= sc.dot(sc.array([[cosphi,sinphi],[-sinphi,cosphi]]),sc.array([pmra,pmdec]))
    return (out[0], out[1])

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
    """
    if sc.array(pmll).shape == ():
        ra,dec = lb_to_radec(l,b,degree=degree,epoch=epoch)
        return pmllpmbb_to_pmrapmdec_single(pmll,pmbb,ra,dec,b,degree,epoch)
    else:
        radec = lb_to_radec(l,b,degree=degree,epoch=epoch)
        function= sc.frompyfunc(pmllpmbb_to_pmrapmdec_single,7,2)
        return sc.array(function(pmll,pmbb,radec[:,0],radec[:,1],b,degree,
                                 epoch)).T

def pmllpmbb_to_pmrapmdec_single(pmll,pmbb,ra,dec,b,degree=False,epoch=2000.0):
    """
    NAME:
       pmllpmbb_to_pmrapmdec_single
    PURPOSE:
       rotate proper motions in (l,b) into proper motions in (ra,dec) for
       scalar inputs
    INPUT:
       pmll - proper motion in l (multplied with cos(b)) [mas/yr]
       pmll - proper motion in b [mas/yr]
       ra - right ascension
       dec - declination
       b - Galactic lattitude
       degree - if True, ra and dec are given in degrees (default=False)
       epoch - epoch of ra,dec (right now only 2000.0 and 1950.0 are supported)
    OUTPUT:
       (pmra,pmdec)
    HISTORY:
       2010-04-07 - Written - Bovy (NYU)
    """
    theta,dec_ngp,ra_ngp= get_epoch_angles(epoch)
    if degree:
        sindec_ngp= m.sin(dec_ngp)
        cosdec_ngp= m.cos(dec_ngp)
        sindec= m.sin(dec*_DEGTORAD)
        sinb= m.sin(b*_DEGTORAD)
        cosdec= m.cos(dec*_DEGTORAD)
        cosb= m.cos(b*_DEGTORAD)
        sinrarangp= m.sin(ra*_DEGTORAD-ra_ngp)
    else:
        sindec_ngp= m.sin(dec_ngp)
        cosdec_ngp= m.cos(dec_ngp)
        sindec= m.sin(dec)
        sinb= m.sin(b)
        cosdec= m.cos(dec)
        cosb= m.sin(b)
        sinrarangp= m.sin(ra*-ra_ngp)
    cosphi= (sindec_ngp-sindec*sinb)/cosdec/cosb
    sinphi= sinrarangp*cosdec_ngp/cosb
    out= sc.dot(sc.array([[cosphi,-sinphi],[sinphi,cosphi]]),sc.array([pmll,pmbb]))
    return (out[0], out[1])


def cov_pmrapmdec_to_pmllpmbb(cov_pmradec,ra,dec,degree=False,epoch=2000.0):
    """
    NAME:
       cov_pmrapmdec_to_pmllpmbb
    PURPOSE:
       propagate the proper motions errors through the rotation from (ra,dec)
       to (l,b)
    INPUT:
       covar_pmradec - uncertainty covariance matrix of the proper motion
                      in ra (multplied with cos(dec)) and dec [2,2] or [:,2,2]
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
        sinb= m.sin(b*_DEGTORAD)
        cosdec= m.cos(dec*_DEGTORAD)
        cosb= m.cos(b*_DEGTORAD)
        sinrarangp= m.sin(ra*_DEGTORAD-ra_ngp)
    else:
        sindec_ngp= m.sin(dec_ngp)
        cosdec_ngp= m.cos(dec_ngp)
        sindec= m.sin(dec)
        sinb= m.sin(b)
        cosdec= m.cos(dec)
        cosb= m.sin(b)
        sinrarangp= m.sin(ra*-ra_ngp)
    cosphi= (sindec_ngp-sindec*sinb)/cosdec/cosb
    sinphi= sinrarangp*cosdec_ngp/cosb
    P= sc.array([[cosphi,sinphi],[-sinphi,cosphi]])
    return sc.dot(P,sc.dot(cov_pmradec,P.T))

def cov_dvrpmllbb_to_vxyz(d,e_d,e_vr,pmll,pmbb,cov_pmllbb,l,b,
                          plx=False,degree=False):
    """
    NAME:
       cov_dvrpmllbb_to_vxyz
    PURPOSE:
       propagate distance, radial velocity, and proper motion uncertainties to
       Galactic coordinates
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
       transform XYZ coordinates (wrt Sun) to rectangular Galactocentric 
       coordinates
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
       transform rectangular Galactocentric to XYZ coordinates (wrt Sun)
       coordinates
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
    if X < 0.:
        phi= m.pi-phi
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
       transform XYZ coordinates (wrt Sun) to cylindrical Galactocentric 
       coordinates
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
       transform cylindrical Galactocentric coordinates to XYZ coordinates 
       (wrt Sun)
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
       transform XYZ coordinates (wrt Sun) to rectangular Galactocentric 
       coordinates for velocities
    INPUT:
       vx - U
       vy - V
       vz - W
       vsun - velocity of the sun ndarray[3]
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
       transform XYZ coordinates (wrt Sun) to cylindrical Galactocentric 
       coordinates for velocities
    INPUT:
       vx - U
       vy - V
       vz - W
       X - X
       Y - Y
       Z - Z
       vsun - velocity of the sun ndarray[3]
       galcen - if True, then X,Y,Z are in cylindrical 
                Galactocentric coordinates
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
       transform rectangular Galactocentric coordinates to XYZ coordinates 
       (wrt Sun) for velocities
    INPUT:
       vXg - Galactocentric x-velocity
       vYg - Galactocentric y-velocity
       vZg - Galactocentric z-velocity
       vsun - velocity of the sun ndarray[3]
    OUTPUT:
       [:,3]= vx, vy, vz
    HISTORY:
       2011-02-24 - Written - Bovy (NYU)
    """
    try:
        return sc.array([-vXg+vsun[0],vYg-vsun[1],vZg-vsun[2]])
    except ValueError: #annoying bug for one-d, make sure they are arrays
        return sc.array([-sc.array([vXg]).flatten()[0]+vsun[0],
                          sc.array([vYg]).flatten()[0]-vsun[1],
                          sc.array([vZg]).flatten()[0]-vsun[2]])      

def galcencyl_to_vxvyvz(vR,vT,vZ,phi,vsun=[0.,1.,0.]):
    """
    NAME:
       galcencyl_to_vxvyvz
    PURPOSE:
       transform cylindrical Galactocentric coordinates to XYZ (wrt Sun)
       coordinates for velocities
    INPUT:
       vR - Galactocentric radial velocity
       vT - Galactocentric tangential velocity
       vZ - Galactocentric vertical velocity
       phi - Galactocentric azimuth
       vsun - velocity of the sun ndarray[3]
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

def get_epoch_angles(epoch=2000.0):
    """
    NAME:
       get_epoch_angles
    PURPOSE:
       get the angles relevant for the transformation from ra, dec to l,b
       for the given epoch
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
        print "Only epochs 1950 and 2000 are supported"
        print "Returning..."
        return -1
    return (theta,dec_ngp,ra_ngp)
