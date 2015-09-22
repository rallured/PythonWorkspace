import numpy as np
import matplotlib.pyplot as plt
import PytranTrace as tran
import zernikemod,pdb,time
##from utilities.plotting import nanmean
import traces.conicsolve as con
import reconstruct

global x,y,z,l,m,n,ux,uy,uz

def transform(tx,ty,tz,rx,ry,rz):
    """Coordinate transformation. translations are done first,
    then Rx,Ry,Rz
    """
    global x,y,z,l,m,n,ux,uy,uz
    tran.transform(x,y,z,l,m,n,ux,uy,uz,-tx,-ty,-tz,-rx,-ry,-rz)
    return

def itransform(tx,ty,tz,rx,ry,rz):
    """Inverse of coordinate transformations. -rz,-ry,-rx then
    translations.
    """
    global x,y,z,l,m,n,ux,uy,uz
    tran.itransform(x,y,z,l,m,n,ux,uy,uz,-tx,-ty,-tz,-rx,-ry,-rz)
    return

def reflect(ind=None):
    """Reflect rays based on surface normal
    """
    global l,m,n,ux,uy,uz
    if ind is not None:
        tl,tm,tn,tux,tuy,tuz = l[ind],m[ind],n[ind],ux[ind],uy[ind],uz[ind]
        tran.reflect(tl,tm,tn,tux,tuy,tuz)
        l[ind],m[ind],n[ind],ux[ind],uy[ind],uz[ind] = tl,tm,tn,tux,tuy,tuz
    else:
        tran.reflect(l,m,n,ux,uy,uz)
    return

def flat(ind=None):
    """Trace rays to the XY plane
    """
    global x,y,z,l,m,n,ux,uy,uz
    if ind is not None:
        tx,ty,tz,tl,tm,tn,tux,tuy,tuz = x[ind],y[ind],z[ind],\
                                        l[ind],m[ind],n[ind],\
                                        ux[ind],uy[ind],uz[ind]
        tran.flat(tx,ty,tz,tl,tm,tn,tux,tuy,tuz)
        x[ind],y[ind],z[ind],\
        l[ind],m[ind],n[ind],\
        ux[ind],uy[ind],uz[ind] = tx,ty,tz,tl,tm,tn,tux,tuy,tuz
    else:
        tran.flat(x,y,z,l,m,n,ux,uy,uz)
    return

def refract(n1,n2):
    """Refract rays based on surface normal
    and ray direction cosines from index n1
    into index n2
    """
    global l,m,n,ux,uy,uz
    tran.refract(l,m,n,ux,uy,uz,n1,n2)
    return

#Wrapper for Zernike surface
#Coordinates are usual arctan2(y,x)
def zernsurf(coeff,rad,rorder=None,aorder=None):
    global x,y,z,l,m,n,ux,uy,uz
    if rorder is None:
        rorder,aorder = zernikemod.zmodes(np.size(coeff))
    tran.tracezern(x,y,z,l,m,n,ux,uy,uz,coeff,\
                   np.array(rorder),np.array(aorder),rad)
    rho = np.sqrt(x**2+y**2)
    ind = np.where(rho<=rad)
    vignette(ind=ind)
    return

#Wrapper for Zernike surface with 2 Zernike sets and one with
#arbitrary rotation angle
#Coordinates are usual arctan2(y,x)
def zernsurfrot(coeff1,coeff2,rad,rot,\
                rorder1=None,aorder1=None,rorder2=None,aorder2=None):
    global x,y,z,l,m,n,ux,uy,uz
    if rorder1 is None:
        rorder1,aorder1 = zernikemod.zmodes(np.size(coeff1))
    if rorder2 is None:
        rorder2,aorder2 = zernikemod.zmodes(np.size(coeff2))
    tran.tracezernrot(x,y,z,l,m,n,ux,uy,uz,coeff1,np.array(rorder1),\
                      np.array(aorder1),coeff2,np.array(rorder2),\
                      np.array(aorder2),rad,rot)
    rho = np.sqrt(x**2+y**2)
    ind = np.where(rho<=rad)
    vignette(ind=ind)

#Wrapper for spherical surface
def sphere(rad):
    global x,y,z,l,m,n,ux,uy,uz
    tran.tracesphere(x,y,z,l,m,n,ux,uy,uz,rad)
    return

def conic(R,K):
    """Wrapper for conic surface with radius of curvature R
    and conic constant K
    """
    global x,y,z,l,m,n,ux,uy,uz
    tran.conic(x,y,z,l,m,n,ux,uy,uz,R,K)
    return

#Wrapper for cylindrical surface
def cyl(rad):
    global x,y,z,l,m,n,ux,uy,uz
    tran.tracecyl(x,y,z,l,m,n,ux,uy,uz,rad)

#Wrapper for cylindrical conics
def cylconic(rad,k):
    global x,y,z,l,m,n,ux,uy,uz
    tran.cylconic(x,y,z,l,m,n,ux,uy,uz,rad,k)

#Wrapper for Wolter primary surface - no vignetting
def wolterprimary(r0,z0):
    global x,y,z,l,m,n,ux,uy,uz
    tran.wolterprimary(x,y,z,l,m,n,ux,uy,uz,r0,z0)
    return

#Wrapper for Wolter secondary surface - no vignetting
def woltersecondary(r0,z0):
    global x,y,z,l,m,n,ux,uy,uz
    tran.woltersecondary(x,y,z,l,m,n,ux,uy,uz,r0,z0)
    return

#Wrapper for Wolter primary surface -
#place at surface tangent point
#+z is surface normal
#+y is toward sky
#+x is azimuthal direction
def wolterprimtan(r0,z0):
    global x,y,z,l,m,n,ux,uy,uz
    #Compute Wolter parameters
    alpha,p,d,e = con.woltparam(r0,z0)
    transform(0,0,0,-np.pi/2-alpha,0,0)
    #Go to Wolter focus minus gap and half mirror length
    transform(0,con.primrad(z0+75.,r0,z0),-z0-75.,0,0,0)
    #Place Wolter surface
    wolterprimary(r0,z0)
    #Go back to original coordinate system
    transform(0,-con.primrad(z0+75.,r0,z0),z0+75.,0,0,0)
    transform(0,0,0,np.pi/2+alpha,0,0)
    return

#Wrapper for Wolter primary surface with sinusoid
def woltersine(r0,z0,amp,freq):
    global x,y,z,l,m,n,ux,uy,uz
    tran.woltersine(x,y,z,l,m,n,ux,uy,uz,r0,z0,amp,freq)
    return

#Wrapper for Wolter sinusoidal surface -
#place at surface tangent point
#+z is surface normal
#+y is toward sky
#+x is azimuthal direction
def woltersinetan(r0,z0,amp,freq):
    global x,y,z,l,m,n,ux,uy,uz
    #Compute Wolter parameters
    alpha,p,d,e = con.woltparam(r0,z0)
    transform(0,0,0,-np.pi/2-alpha,0,0)
    #Go to Wolter focus minus gap and half mirror length
    transform(0,con.primrad(z0+75.,r0,z0),-z0-75.,0,0,0)
    #Place Wolter surface
    woltersine(r0,z0,amp,freq)
    #Go back to original coordinate system
    transform(0,-con.primrad(z0+75.,r0,z0),z0+75.,0,0,0)
    transform(0,0,0,np.pi/2+alpha,0,0)
    return

##Wrapper for L-L secondary surface
##Placed at focus
def secondaryLL(r0,z0,zmax,zmin,dphi,coeff,axial,az):
    global x,y,z,l,m,n,ux,uy,uz
    tran.woltersecll(x,y,z,l,m,n,ux,uy,uz,r0,z0,zmax,zmin,dphi,coeff,axial,az)
    vignette()
    return

##Wrapper for L-L primary surface
##Placed at focus
def primaryLL(r0,z0,zmax,zmin,dphi,coeff,axial,az):
    global x,y,z,l,m,n,ux,uy,uz
    tran.wolterprimll(x,y,z,l,m,n,ux,uy,uz,r0,z0,zmax,zmin,dphi,coeff,axial,az)
    vignette()
    return

#Wrapper for Wolter primary surface -
#place at surface tangent point
#+z is surface normal
#+y is toward sky
#+x is azimuthal direction
def primaryLLtan(r0,z0,zmax,zmin,dphi,coeff,axial,az):
    global x,y,z,l,m,n,ux,uy,uz
    #Compute Wolter parameters
    alpha,p,d,e = con.woltparam(r0,z0)
    transform(0,0,0,-np.pi/2-alpha,0,0)
    #Go to Wolter focus minus gap and half mirror length
    transform(0,con.primrad(z0+75.,r0,z0),-z0-75.,0,0,0)
    #Place Wolter surface
    transform(0,0,0,0,0,-np.pi/2)
    primaryLL(r0,z0,zmax,zmin,dphi,coeff,axial,az)
    transform(0,0,0,0,0,np.pi/2)
    #Go back to original coordinate system
    transform(0,-con.primrad(z0+75.,r0,z0),z0+75.,0,0,0)
    transform(0,0,0,np.pi/2+alpha,0,0)
    return

#Wrapper for Wolter-Schwarzschild Primary
def wsPrimary(r0,z0,psi):
    """Trace a W-S primary surface
    Fortran function computes Chase parameters for an equivalent W-I
    betas, f, g, and k computed from alpha and z0
    """
    global x,y,z,l,m,n,ux,uy,uz
    a,p,d,e = con.woltparam(r0,z0)
    tran.wsprimary(x,y,z,l,m,n,ux,uy,uz,a,z0,psi)
    return

#Wrapper for Wolter-Schwarzschild Secondary
def wsSecondary(r0,z0,psi):
    """Trace a W-S secondary surface
    Fortran function computes Chase parameters for an equivalent W-I
    betas, f, g, and k computed from alpha and z0
    """
    global x,y,z,l,m,n,ux,uy,uz
    a,p,d,e = con.woltparam(r0,z0)
    tran.wssecondary(x,y,z,l,m,n,ux,uy,uz,a,z0,psi)
    return

#Wrapper for SPO surface
def spoCone(R0,tg,ind=None):
    """Trace rays to an SPO cone with intersection radius
    R0 and slope angle tg.
    XY plane should be at SPO intersection plane
    """
    global x,y,z,l,m,n,ux,uy,uz
    if ind is not None:
        tx,ty,tz,tl,tm,tn,tux,tuy,tuz = x[ind],y[ind],z[ind],\
                                        l[ind],m[ind],n[ind],\
                                        ux[ind],uy[ind],uz[ind]
        tran.spocone(tx,ty,tz,tl,tm,tn,tux,tuy,tuz,R0,tg)
        x[ind],y[ind],z[ind],\
        l[ind],m[ind],n[ind],\
        ux[ind],uy[ind],uz[ind] = tx,ty,tz,tl,tm,tn,tux,tuy,tuz
    else:
        tran.spocone(x,y,z,l,m,n,ux,uy,uz,R0,tg)
    return

#Wrapper for SPO as function of radius and focal length
def spoPrimary(R0,F,d=.605,ind=None):
    """Trace rays to an SPO primary with intersection radius
    R0 and focal length F.
    XY plane should be at SPO intersection plane
    """
    #Convert F to tg
    tg = .25*np.arctan((R0+d/2)/F)
    #Call SPO wrapper
    spoCone(R0,tg,ind=ind)
    return

#Wrapper for SPO as function of radius and focal length
def spoSecondary(R0,F,d=.605,ind=None):
    """Trace rays to an SPO secondary with intersection radius
    R0 and focal length F.
    XY plane should be at SPO intersection plane
    """
    #Convert F to tg
    tg = .75*np.arctan((R0+d/2)/F)
    #Call SPO wrapper
    spoCone(R0,tg,ind=ind)
    return

#Wrapper for radial grating
def radgrat(hubdist,dpermm,order,wave,ind=None):
    """Infinite radial grating. Assumes grating in x,y plane
    with grooves converging at hubdist in positive y direction
    dpermm is nm/mm
    wave is in nm
    """
    global x,y,l,m,n
    if ind is not None:
        tx,ty,tl,tm,tn = x[ind],y[ind],l[ind],m[ind],n[ind]
        tran.radgrat(tx,ty,tl,tm,tn,hubdist,dpermm,order,wave)
        x[ind],y[ind],l[ind],m[ind],n[ind] = tx,ty,tl,tm,tn
    else:
        tran.radgrat(x,y,l,m,n,hubdist,dpermm,order,wave)
    return

def grat(d,order,wave):
    """Linear grating with groove direction in +y
    Evanescence results in position vector set to zero
    """
    tran.grat(x,y,l,m,n,d,order,wave)
    return

#Remove vignetted rays from memory
#ind is array of "good" indices, all others are removed
def vignette(ind=None):
    global x,y,z,l,m,n,ux,uy,uz
    if ind==None:
        mag = l**2+m**2+n**2
        ind = np.where(mag>.1) #Automatic vignetting
                            #requires position vector set to 0.
    x = x[ind]
    y = y[ind]
    z = z[ind]
    l = l[ind]
    m = m[ind]
    n = n[ind]
    ux = ux[ind]
    uy = uy[ind]
    uz = uz[ind]
    return

#Trace lens, first surface center is coincident with xy plane
#thickness extends in positive z direction
#Rays should be traced to this plane before calling lens
#to ensure rays hit correct spherical intersection
def lens(r1,r2,thick,d,nl,reverse=False):
    global x,y,z,l,m,n,ux,uy,uz
    if reverse is True:
        r1,r2 = -r2,-r1
    #Trace to first surface
    if r1 != 0:
        transform(0.,0.,r1,0.,0.,0.) #Go to first
                                     #center of curvature, r1>0->convex
        sphere(abs(r1)) #Trace to first spherical surface
    rho = np.sqrt(x**2 + y**2)/(d/2)
    ind = np.where(rho<=1.) #Remove rays with radius beyond that of lens
    vignette(ind=ind)
    refract(1.,nl) #Refract into lens

    
    transform(0.,0.,-r1+thick,0.,0.,0.) #Go to center of second surface
    flat()
    if r2 != 0:
        transform(0.,0.,r2,0,0,0) #Go to center of curvature
        sphere(abs(r2)) #Trace to second spherical surface
    rho = np.sqrt(x**2 + y**2)/(d/2)
    ind = np.where(rho<=1.) #Remove rays with radius beyond that of lens
    vignette(ind=ind)
    refract(nl,1.) #Refract out of lens
    transform(0.,0.,-r2,0,0,0) #Transform to xy plane tangent to
                              #center of second surface

    #Rays are left at second surface, pointing out of lens in correct direction
    return

#Cylindrical lens, same principle as with standard lens
def cyllens(r1,r2,thick,width,height,nl,reverse=False):
    global x,y,z,l,m,n,ux,uy,uz
    if reverse is True:
        r1,r2 = -r2,-r1
    #Trace to first surface
    if r1 != 0:
        transform(0.,0,r1,0.,0.,0.) #Go to first
                                     #center of curvature, r1>0->convex
        cyl(abs(r1)) #Trace to first spherical surface
    ind = logical_and(x<width/2,y<height/2) #Remove rays outside of cylinder
    vignette(ind=ind)
    refract(1.,nl) #Refract into lens

    
    transform(0.,0,-r1+thick,0.,0.,0.) #Go to center of second surface
    flat()
    if r2 != 0:
        transform(0.,0,r2,0,0,0) #Go to center of curvature
        cyl(abs(r2)) #Trace to second spherical surface
    ind = logical_and(x<width/2,y<height/2) #Remove rays outside of cylinder
    vignette(ind=ind)
    refract(nl,1.) #Refract out of lens
    transform(0.,0,-r2,0,0,0) #Transform to xy plane tangent to
                              #center of second surface

    #Rays are left at second surface, pointing out of lens in correct direction
    return


###### SOURCES #######

#Define point source with angular divergence
#Points in +z direction
def pointsource(ang,num):
    global x,y,z,l,m,n,ux,uy,uz
    #Radial direction cosine magnitude
    rho = np.sqrt(np.random.rand(num))*np.sin(ang)
    theta = np.random.rand(num)*2*np.pi
    l = rho*np.cos(theta)
    m = rho*np.sin(theta)
    n = np.sqrt(1.-l**2-m**2)
    x = np.repeat(0.,num)
    y = np.repeat(0.,num)
    z = np.repeat(0.,num)
    ux = np.repeat(0.,num)
    uy = np.repeat(0.,num)
    uz = np.repeat(0.,num)
    return

#Define uniform, circular beam of radius rad, pointing in +z direction
def circularbeam(rad,num):
    global x,y,z,l,m,n,ux,uy,uz
    rho = np.sqrt(np.random.rand(num))*rad
    theta = np.random.rand(num)*2*np.pi
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(0.,num)
    l = np.repeat(0.,num)
    m = np.repeat(0.,num)
    n = np.repeat(1.,num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    return

#Define annulus of rays pointing in +z direction
def annulus(rin,rout,num):
    global x,y,z,l,m,n,ux,uy,uz
    rho = np.sqrt(rin**2+np.random.rand(num)*(rout**2-rin**2))
    theta = np.random.rand(num)*2*np.pi
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(0.,num)
    l = np.repeat(0.,num)
    m = np.repeat(0.,num)
    n = np.repeat(1.,num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    return

def subannulus(rin,rout,dphi,num):
    """Create a subapertured annulus source in +z direction
    Annulus is centered about theta=0 which points to +x
    """
    global x,y,z,l,m,n,ux,uy,uz
    rho = np.sqrt(rin**2+np.random.rand(num)*(rout**2-rin**2))
    theta = np.random.rand(num)*dphi - dphi/2.
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(0.,num)
    l = np.repeat(0.,num)
    m = np.repeat(0.,num)
    n = np.repeat(1.,num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    return

#Converging beam
#Place at nominal focus
#Input z position, inner and outer radius,
#min and max theta
def convergingbeam(zset,rin,rout,tmin,tmax,num,lscat):
    global x,y,z,l,m,n,ux,uy,uz
    rho = sqrt(rin**2+random.rand(num)*(rout**2-rin**2))
    theta = tmin + random.rand(num)*(tmax-tmin)
    x = rho*cos(theta)
    y = rho*sin(theta)
    z = np.repeat(zset,num)
    lscat = lscat * tan((random.rand(num) - .5)*np.pi)
    lscat = lscat/60**2 * np.pi/180.
    n = -cos(arctan(rho/zset)+lscat)
    l = -sqrt(1-n**2)*cos(theta)
    m = -sqrt(1-n**2)*sin(theta)
    ux = np.repeat(0.,num)
    uy = np.repeat(0.,num)
    uz = np.repeat(0.,num)
    return

#Converging beam
#Place at nominal focus
#Input z position, inner and outer radius,
#min and max theta
def convergingbeam2(zset,xmin,xmax,ymin,ymax,num,lscat):
    global x,y,z,l,m,n,ux,uy,uz
    x = xmin + random.rand(num)*(xmax-xmin)
    y = ymin + random.rand(num)*(ymax-ymin)
    rho = sqrt(x**2+y**2)
    theta = arctan2(y,x)
    z = np.repeat(zset,num)
    lscat = lscat * tan((random.rand(num) - .5)*np.pi)
    lscat = lscat/60**2 * np.pi/180.
    n = -cos(arctan(rho/zset)+lscat)
    l = -sqrt(1-n**2)*cos(theta)
    m = -sqrt(1-n**2)*sin(theta)
    ux = np.repeat(0.,num)
    uy = np.repeat(0.,num)
    uz = np.repeat(0.,num)
    return

#Rectangular beam pointing in +z direction
def rectbeam(xhalfwidth,yhalfwidth,num):
    global x,y,z,l,m,n,ux,uy,uz
    x = (random.rand(num)-.5)*2*xhalfwidth
    y = (random.rand(num)-.5)*2*yhalfwidth
    z = np.repeat(0.,num)
    n = np.repeat(1.,num)
    l = np.copy(z)
    m = np.copy(z)
    ux = np.repeat(0.,num)
    uy = np.repeat(0.,num)
    uz = np.repeat(0.,num)
    return
    
####### Lenses ##########
def LJ1653L2(reverse=False):
    cyllens(103.36,0,4.09,30.,60.,1.51501,reverse=reverse)
    return

def LJ1629L2(reverse=False):
    cyllens(77.52,0,4.46,30.,60.,1.51501,reverse=reverse)
    return

def AC254_400_A(reverse=False):
    if reverse is False:
        lens(738.5,181.55,2.,12.7*2,1.64363)
        lens(181.55,-219.8,4.,12.7*2,1.51501)
    else:
        lens(181.55,-219.8,4.,12.7*2,1.51501,reverse=True)
        lens(738.5,181.55,2.,12.7*2,1.64363,reverse=True)
    return

def AC508_200_A(reverse=False):
    if reverse is False:
        lens(109.86,-93.110,8.5,50.8,1.51501)
        lens(-93.110,-376.25,2.,50.8,1.64363)
    else:
        lens(-93.110,-376.25,2.,50.8,1.64363,reverse=True)
        lens(109.86,-93.110,8.5,50.8,1.51501,reverse=True)
    return

def cylNull(reverse=False):
    if reverse is False:
        transform(0,0,0,np.pi/2,0,0)
        cylconic(.007626,-.575)
        refract(1.,1.51501)
        transform(0,0,0,-np.pi/2,0,0)
        transform(0,0,50,0,0,0)
        flat()
        refract(1.51501,1.)
    else:
        refract(1.,1.51501)
        transform(0,0,50,0,0,0)
        transform(0,0,0,np.pi/2,0,0)
        cylconic(-.007626,-.575)
        refract(1.51501,1.)
        transform(0,0,0,-np.pi/2,0,0)
    return


#######  ANALYSES #########
def centroid(weights=None):
    """Compute the centroid of the rays in the xy plane
    """
    cx = np.average(x,weights=weights)
    cy = np.average(y,weights=weights)
    return cx,cy

def rmsCentroid(weights=None):
    """Compute the RMS of rays from centroid in xy plane
    """
    cx,cy = centroid(weights=weights)
    rho = (x-cx)**2 + (y-cy)**2
    return np.sqrt(np.average(rho,weights=weights))

def hpd(weights=None):
    """Compute HPD by taking median of radii from centroid"""
    cx,cy = centroid(weights=weights)
    rho = np.sqrt((x-cx)**2 + (y-cy)**2)
    if weights is not None:
        ind = np.argsort(rho)
        weights = weights[ind]
        rho = rho[ind]
        cdf = np.cumsum(weights)
        hpd = rho[np.argmin(np.abs(cdf-.5))] * 2.
    else:
        hpd = np.median(rho)*2.
    return hpd

def findimageplane(zscan,num,weights=None):
    global x,y,z,l,m,n,ux,uy,uz
    rms = []
    zsteps = np.linspace(-zscan,zscan,num)
    for znow in np.linspace(-zscan,zscan,num):
        #Transform to offset
        transform(0,0,znow,0,0,0)
        #Trace rays to new plane
        flat()
        rms.append(rmsCentroid(weights=weights))
        #Return to nominal plane
        transform(0,0,-znow,0,0,0)
    flat()

##    plt.clf()
##    plt.plot(zsteps,rms)

    return zsteps[np.where(rms==np.min(rms))]

def findlineplane(zscan,num):
    global x,y,z,l,m,n,ux,uy,uz
    rms = []
    zsteps = linspace(-zscan,zscan,num)
    for znow in linspace(-zscan,zscan,num):
        #Transform to offset
        transform(0,0,znow,0,0,0)
        #Trace rays to new plane
        flat()
        #Determine centroid
        cx = nanmean(x)
        #Compute RMS image size
        rho2 = (x-cx)**2
        rms.append(sqrt(mean(rho2)))
        #Return to nominal plane
        transform(0,0,-znow,0,0,0)
    flat()

    clf()
    plot(zsteps,rms)

    return zsteps[where(rms==min(rms))]

def findimageplane2(zscan,num):
    global x,y,z,l,m,n,ux,uy,uz
    hew = []
    zsteps = linspace(-zscan,zscan,num)
    for znow in linspace(-zscan,zscan,num):
        #Transform to offset
        transform(0,0,znow,0,0,0)
        #Trace rays to new plane
        flat()
        #Determine centroid
        cx = nanmean(x)
        cy = nanmean(y)
        #Compute RMS image size
        rho2 = sqrt((x-cx)**2+(y-cy)**2)
        hew.append(median(rho2)*2.)
        #Return to nominal plane
        transform(0,0,-znow,0,0,0)
    flat()

    clf()
    plot(zsteps,hew)

    return zsteps[where(hew==min(hew))]

def wsPrimRad(z,psi,r0,z0):
    """Return the radius of a WS primary as a function of axial coordinate
    This is computed numerically by tracing a single ray in plane
    orthogonal to optical axis
    """
    #Set up source pointing toward +z
    pointsource(0.,1)
    transform(0,0,0,0,-np.pi/2,0) #Point ray to +x
    transform(-r0,0,-z,0,0,0) #Go to proper axial location

    #Trace to WS primary
    wsPrimary(r0,z0,psi)

    return x[0]

def wsSecRad(z,psi,r0,z0):
    """Return the radius of a WS primary as a function of axial coordinate
    This is computed numerically by tracing a single ray in plane
    orthogonal to optical axis
    """
    #Set up source pointing toward +z
    pointsource(0.,1)
    transform(0,0,0,0,-np.pi/2,0) #Point ray to +x
    transform(-r0,0,-z,0,0,0) #Go to proper axial location

    #Trace to WS primary
    wsSecondary(r0,z0,psi)

    return x[0]

#Wrapper to create array of zernike values
def zerntest(n,m):
    #Plot radialpoly vs. fastradpoly for rho=0 to 1
    rho = linspace(0.,1.,10)
    vec1 = np.copy(rho)
    vec2 = np.copy(rho)
    for i in range(size(vec1)):
        vec1[i] = tran.radialpoly(rho[i],n,m)
        vec2[i] = tran.fastradpoly(rho[i],n,m)
    plot(vec1,vec2,'.')
    return

#Wrapper for reconstructing referenced wavefront
def referencedWavefront(xang,yang,phase,xang2,yang2,phase2):
    phaseinf = np.copy(phase)
    phaseinf[:,:] = 0.
    ind = where(logical_or(phase==100,phase2==100))
    phaseinf[ind] = 100
    xanginf = np.copy(xang)
    xanginf = xang2-xang
    yanginf = np.copy(yang)
    yanginf = yang2-yang
    xanginf[ind] = 100
    yanginf[ind] = 100

    #Reconstruct influence wavefront
    influence = reconstruct.reconstruct(xanginf,yanginf,1.e-12,phaseinf)

    #Make invalid np.pixels NaNs
    ind = where(influence==100.)
    influence[ind] = NaN

    return influence
