import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import pdb
from utilities.imaging.fitting import circle,circleMerit
from utilities.imaging.analysis import ptov,rms
import utilities.imaging.man as man

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
    cy = np.average(y,weights=weights)
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

def grazeAngle(rays,flat=False):
    """Find the graze angle of the rays with the current
    surface normal."""
    return np.arcsin(rays[4]*rays[7] +\
                     rays[5]*rays[8] +\
                     rays[6]*rays[9])

def interpolateVec(rays,I,Nx,Ny,xr=None,yr=None,method='linear'):
    """
    Interpolate a ray vector onto a 2D grid based on the X and Y
    positions of the rays. Assume that the rays randomly fill a
    contiguous region of space. Automatically set up the grid
    using the min/max X/Y positions. User inputs the number of
    points in the grid as Nx,Ny, where total points is Nx*Ny.
    User can also indicate the interpolation method used by
    griddata.
    I indicates which vector to interpolate
    """
    #Unpack needed vectors
    x,y = rays[1:3]
    interpVec = rays[I]
    #Set up new grid
    if xr is None:
        xr=[x.min(),x.max()]
        yr=[y.min(),y.max()]
    gridx,gridy = np.meshgrid(np.linspace(xr[0],xr[1],Nx),\
                              np.linspace(yr[0],yr[1],Ny))
    dx = np.diff(gridx)[0][0]
    dy = np.diff(np.transpose(gridy))[0][0]
    #Call griddata for interpolation
    res = griddata((x,y),interpVec,(gridx,gridy),method=method)
    
    return res,dx,dy

def measurePower(rays,Nx,Ny,method='linear'):
    """Measure the radius of curvature in X and Y
    axes of OPD
    Assumes you have steered the beam to get rid of
    any tilts
    Sign of power equals transform in z to go to focus"""
    #Get l,m
    l,dx,dy = interpolateVec(rays,4,Nx,Ny,method=method)
    m = interpolateVec(rays,5,Nx,Ny,method=method)[0]
    #Get slices
    xsl = man.stripnans(l[Ny/2])
    ysl = man.stripnans(m[:,Nx/2])
    #Estimate gradients
    xpow = 1/np.gradient(xsl,dx)[Nx/2]
    ypow = 1/np.gradient(ysl,dy)[Ny/2]
    return -xpow,-ypow
    

def compareOPDandSlopes(rays,Nx,Ny,method='linear'):
    """
    Function to compare OPD gradient to ray slopes quickly
    """
    #Get opd
    opd,dx = interpolateVec(rays,0,Nx,Ny,method=method)
    #Get gradient
    gradx,grady = np.gradient(opd,dx)
    #Get slopes
    l,dx = interpolateVec(rays,4,Nx,Ny,method=method)
    m,dx = interpolateVec(rays,5,Nx,Ny,method=method)
    #Make plots
    fig = plt.figure()
    
    fig.add_subplot(321)
    plt.imshow(grady)
    plt.title('Gradx')
    plt.colorbar()
    
    fig.add_subplot(322)
    plt.imshow(gradx)
    plt.title('Grady')
    plt.colorbar()

    fig.add_subplot(323)
    plt.imshow(l)
    plt.title('l')
    plt.colorbar()

    fig.add_subplot(324)
    plt.imshow(m)
    plt.title('m')
    plt.colorbar()

    fig.add_subplot(325)
    plt.imshow(l-grady)
    plt.title('l-gradx')
    plt.colorbar()

    fig.add_subplot(326)
    plt.imshow(m-gradx)
    plt.title('m-grady')
    plt.colorbar()

    print 'X Diff:' + str(np.sqrt(np.nanmean((grady-l)**2)))
    print 'Y Diff:' + str(np.sqrt(np.nanmean((gradx-m)**2)))
    
    return
