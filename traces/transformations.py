import numpy as np
import transformationsf as tran


def transform(rays,tx,ty,tz,rx,ry,rz,ind=None):
    """Coordinate transformation. translations are done first,
    then Rx,Ry,Rz
    """
    x,y,z,l,m,n,ux,uy,uz = rays[1:]
    if ind is not None:
        tx,ty,tz,tl,tm,tn,tux,tuy,tuz = x[ind],y[ind],z[ind],\
                                        l[ind],m[ind],n[ind],\
                                        ux[ind],uy[ind],uz[ind]
        tran.transform(tx,ty,tz,tl,tm,tn,tux,tuy,tuz,-tx,-ty,-tz,-rx,-ry,-rz)
        x[ind],y[ind],z[ind],\
        l[ind],m[ind],n[ind],\
        ux[ind],uy[ind],uz[ind] = tx,ty,tz,tl,tm,tn,tux,tuy,tuz
    else:
        tran.transform(x,y,z,l,m,n,ux,uy,uz,-tx,-ty,-tz,-rx,-ry,-rz)
    return


def itransform(rays,tx,ty,tz,rx,ry,rz):
    """Inverse of coordinate transformations. -rz,-ry,-rx then
    translations.
    """
    x,y,z,l,m,n,ux,uy,uz = rays[1:]
    tran.itransform(x,y,z,l,m,n,ux,uy,uz,-tx,-ty,-tz,-rx,-ry,-rz)
    return

def reflect(rays,ind=None):
    """Reflect rays based on surface normal
    """
    l,m,n,ux,uy,uz = rays[4:]
    if ind is not None:
        tl,tm,tn,tux,tuy,tuz = l[ind],m[ind],n[ind],ux[ind],uy[ind],uz[ind]
        tran.reflect(tl,tm,tn,tux,tuy,tuz)
        l[ind],m[ind],n[ind],ux[ind],uy[ind],uz[ind] = tl,tm,tn,tux,tuy,tuz
    else:
        tran.reflect(l,m,n,ux,uy,uz)
    return

def refract(rays,n1,n2):
    """Refract rays based on surface normal
    and ray direction cosines from index n1
    into index n2
    """
    l,m,n,ux,uy,uz = rays[4:]
    tran.refract(l,m,n,ux,uy,uz,n1,n2)
    return

def radgrat(rays,hubdist,dpermm,order,wave,ind=None):
    """Infinite radial grating. Assumes grating in x,y plane
    with grooves converging at hubdist in positive y direction
    dpermm is nm/mm
    wave is in nm
    """
    x,y,z,l,m,n = rays[1:7]
    if ind is not None:
        tx,ty,tl,tm,tn = x[ind],y[ind],l[ind],m[ind],n[ind]
        tran.radgrat(tx,ty,tl,tm,tn,hubdist,dpermm,order,wave)
        x[ind],y[ind],l[ind],m[ind],n[ind] = tx,ty,tl,tm,tn
    else:
        tran.radgrat(x,y,l,m,n,hubdist,dpermm,order,wave)
    return

def grat(rays,d,order,wave):
    """Linear grating with groove direction in +y
    Evanescence results in position vector set to zero
    """
    x,y,z,l,m,n = rays[1:7]
    tran.grat(x,y,l,m,n,d,order,wave)
    return

def vignette(rays,ind=None):
    """Remove vignetted rays from memory
    ind is array of "good" indices, all others are removed
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if ind==None:
        mag = l**2+m**2+n**2
        ind = np.where(mag>.1) #Automatic vignetting
                            #requires position vector set to 0.
    
    return [rays[i][ind] for i in range(10)]
