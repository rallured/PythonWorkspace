import numpy as np

def centroid(rays,weights=None):
    """Compute the centroid of the rays in the xy plane
    """
    x,y = rays[1:3]
    cx = np.average(x,weights=weights)
    cy = np.average(y,weights=weights)
    return cx,cy

def rmsCentroid(rays,weights=None):
    """Compute the RMS of rays from centroid in xy plane
    """
    x,y = rays[1:3]
    cx,cy = centroid(rays,weights=weights)
    rho = (x-cx)**2 + (y-cy)**2
    return np.sqrt(np.average(rho,weights=weights))

def rmsX(rays,weights=None):
    """RMS from centroid in the X direction"""
    x = rays[1]
    cx = np.average(x,weights=weights)
    rmsx = np.sqrt(np.average((x-cx)**2,weights=weights))
    return rmsx

def rmsY(rays,weights=None):
    """RMS from centroid in the Y direction"""
    y = rays[2]
    cy = np.average(y,weights=weights)
    rmsy = np.sqrt(np.average((y-cy)**2,weights=weights))
    return rmsy

def hpd(rays,weights=None):
    """Compute HPD by taking median of radii from centroid"""
    x,y = rays[1:3]
    cx,cy = centroid(rays,weights=weights)
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

def hpdY(rays,weights=None):
    """Compute HPD in y direction by taking median of radii from centroid
    Does rho need to be absolute value???
    """
    y = rays[2]
    cy = np.average(PT.y,weights=weights)
    rho = np.abs(y-cy)
    if weights is not None:
        ind = np.argsort(rho)
        weights = weights[ind]
        rho = rho[ind]
        cdf = np.cumsum(weights)
        hpd = rho[np.argmin(np.abs(cdf-.5))] * 2.
    else:
        hpd = np.median(rho)*2.
    return hpd

def analyticImagePlane(rays,weights=None):
    """Find the image plane using the analytic method from
    Ron Elsner's paper
    """
    x,y,z,l,m,n = rays[1:7]
    bx = np.average(x*l/n,weights=weights)-np.average(x,weights=weights)\
         *np.average(l/n,weights=weights)
    ax = np.average((l/n)**2,weights=weights)\
         -np.average(l/n,weights=weights)**2
    by = np.average(y*m/n,weights=weights)-np.average(y,weights=weights)\
         *np.average(m/n,weights=weights)
    ay = np.average((m/n)**2,weights=weights)\
         -np.average(m/n,weights=weights)**2
    dz = -(bx+by)/(ax+ay)
    
    return dz

def analyticYPlane(rays,weights=None):
    """Find the line plane using analytic method from
    Ron Elsner's paper"""
    x,y,z,l,m,n = rays[1:7]
    by = np.average(y*m/n,weights=weights)-np.average(y,weights=weights)\
         *np.average(m/n,weights=weights)
    ay = np.average((m/n)**2,weights=weights)\
         -np.average(m/n,weights=weights)**2
    dz = -by/ay
    return dz

def analyticXPlane(rays,weights=None):
    """Find the line plane using analytic method from
    Ron Elsner's paper"""
    x,y,z,l,m,n = rays[1:7]
    bx = np.average(x*l/n,weights=weights)-np.average(x,weights=weights)\
         *np.average(l/n,weights=weights)
    ax = np.average((l/n)**2,weights=weights)\
         -np.average(l/n,weights=weights)**2
    dz = -bx/ax
    return dz

def grazeAngle(rays):
    """Find the graze angle of the rays with the current
    surface normal."""
    return np.arcsin(rays[4]*rays[7] +\
                     rays[5]*rays[8] +\
                     rays[6]*rays[9])
