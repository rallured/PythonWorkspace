import numpy as np
import matplotlib.pyplot as plt
import pyfits
from scipy.interpolate import griddata
import pdb
import traces.conicsolve as con
import traces.PyTrace as PT
import pyfits
from utilities.imaging.analysis import rms
import utilities.imaging.fitting as fit

#Global directory variables for problem
datadir = '/home/rallured/data/solve_pzt/'

def convertOABInfluence(filename,Nx,Ny,method='cubic'):
    """Read in Vanessa's CSV file for 1 m RoC AXRO mirror
    Convert the azimuthal cartesian coordinate into an arc length
    Regrid the perturbed coordinates onto a regular arc and axial
    coordinate grid
    """
    d = np.transpose(np.genfromtxt(filename,skip_header=1,delimiter=','))
##    x0,y0,z0 = pyfits.getdata('150519StartingNodes.fits')
##    x = x0 + d[5]
##    y = y0 + d[6]
##    z = z0 + d[7]
    x = d[2] + d[5]
    y = d[3] + d[6]
    z = d[4] + d[7]
    pdb.set_trace()

    theta = np.arctan2(y,-z)
    rPert = 1. - np.sqrt(x**2+z**2) #Positive result is toward OA
    arc = theta * 1.

    #Now have x for cylindrical axis, arc for azimuthal axis,
    #and rPert for the radial perturbation
    #Regrid to a regular grid and you're done
    gy = np.linspace(x.min(),x.max(),Nx+2)
    gx = np.linspace(arc.min(),arc.max(),Ny+2)
    gx,gy = np.meshgrid(gx,gy)
    d = np.transpose(griddata((x,arc),rPert,(gy,gx),method=method))
    d[np.isnan(d)] = 0.

    return d[1:-1,1:-1],gx,gy

def oabDistortion(amp,freq,phase,filename):
    """Similar to the original createDistortion, this introduces
    a sinusoidal ripple to the 1 meter RoC geometry.
    """
    #Create distortion array
    y,x = np.mgrid[0:200,0:200]
    x,y = x*.5,y*.5
    d = amp*np.sin(2*np.pi*freq*y+phase)

    #Save as fits file
    pyfits.writeto(datadir+'distortions/'+filename,d,clobber=True)

    return d

def wolterDistortion(filename,Nx,Ny):
    """This creates the distortion map of a 1 meter RoC
    Wolter I as compared to a 1 meter cylinder
    Should be matched to height/theta grid
    Of influence functions
    """
    #Determine Wolter prescription by finding
    #intersection radius closest to 1 meter
    #at 1125 mm above focal point
    rint = np.linspace(990.,1000.,1000)
    r0 = con.primrad(10125.,rint,10000.)
    rselect = rint[np.argmin(abs(r0-1000.))]
    a,p,d,e = con.woltparam(rselect,10000.)
    
    #Need to get rays onto Wolter I surface
    tspan = .04952303577 #maximum theta from IFs
    zspan = 50.*np.cos(2*a)
    t,z = np.meshgrid(np.linspace(-tspan,tspan,200),\
                      np.linspace(10125.-zspan,10125.+zspan,200))
    #Set up source rays along optical axis
    PT.z = z.flatten()
    PT.l = np.cos(t.flatten())
    PT.m = np.sin(t.flatten())
    PT.n = np.zeros(np.size(PT.z))
    PT.x = PT.l*1000.
    PT.y = PT.m*1000.
    PT.ux = np.zeros(np.size(PT.z))
    PT.uy = np.zeros(np.size(PT.z))
    PT.uz = np.zeros(np.size(PT.z))
    #Trace to Wolter I surface and go to cylinder
    PT.wolterprimary(rselect,10000.)
    PT.transform(con.primrad(10125.,rselect,10000.),0,10125.,0,0,0)
    PT.transform(0,0,0,0,a,0)
    PT.transform(-1000.,0,0,0,0,0)
    pdb.set_trace()

    #Compute radial perturbations, reshape, and return distortion
    rpert = 1000.-np.sqrt(PT.x**2+PT.y**2)
    return rpert.reshape((200,200)),z.reshape((200,200)),t.reshape((200,200))

def determineHFDFC2Stress(fn):
    """Determine required stress*thickness for correctin of
    HFDFDC2 after PZT processing."""
##    '/home/rallured/Dropbox/AXRO/HFDFC/'
##                        '110x110_50x250_200Hz_xyscan_'
##                        'Height_transformed_4in_deltaR_matrix.fits'
    #Load in data and create coordinate grids
    fig = pyfits.getdata(fn)/1e6
    xf,yf = np.meshgrid(np.linspace(-.045,.045,385),\
                        np.linspace(-.045,.045,385))
    #Fill in NaNs with linear interpolation
    xf2,yf2 = xf[~np.isnan(fig)],yf[~np.isnan(fig)]
    fig2 = griddata((xf2,yf2),fig[~np.isnan(fig)],(xf,yf))
    fig2[0,0] = fig2[0,1]
    #Perform SG smoothing
    gf = fit.sgolay2d(fig2,41,3,derivative='col')
    #Load stress data
    uni = pyfits.getdata('/home/rallured/Dropbox/AXRO/HFDFC/Uni.fits')
    xu,yu = np.meshgrid(np.linspace(-.05,.05,200),\
                        np.linspace(-.05,.05,200))
    #Interpolate stress onto figure grid
    uni2 = griddata((xu.flatten(),yu.flatten()),uni.flatten(),\
                       (xf,yf))
    #Compute gradients
    dx = .09/(399) #Grid size in meters
    gf = gf/dx
    gs = np.gradient(uni2,dx)[0]
    #Remove tilt
    gf = gf - np.nanmean(gf)
    gs = gs - np.nanmean(gs)
    #Determine coefficient needed to match axial sags
    fun = lambda c: rms(gf-c*gs)
    coeff = np.linspace(-5,5,1000)
    fom = [fun(c) for c in coeff]
    pdb.set_trace()

def determineHFDFC2StressSurf():
    """Determine required stress*thickness for correctin of
    HFDFDC2 after PZT processing."""
    #Load in data and create coordinate grids
    fig = pyfits.getdata('/home/rallured/Dropbox/AXRO/HFDFC/'
                        '110x110_50x250_200Hz_xyscan_'
                        'Height_transformed_4in_deltaR_matrix.fits')/1e6
    xf,yf = np.meshgrid(np.linspace(-.045,.045,400),\
                        np.linspace(-.045,.045,400))
    #Fill in NaNs with linear interpolation
    xf2,yf2 = xf[~np.isnan(fig)],yf[~np.isnan(fig)]
    fig2 = griddata((xf2,yf2),fig[~np.isnan(fig)],(xf,yf))
    #Perform SG smoothing
    fig2 = fit.sgolay2d(fig2,41,3)
    #Load stress data
    uni = pyfits.getdata('/home/rallured/Dropbox/AXRO/HFDFC/Uni.fits')
    xu,yu = np.meshgrid(np.linspace(-.05,.05,200),\
                        np.linspace(-.05,.05,200))
    #Interpolate stress onto figure grid
    uni2 = griddata((xu.flatten(),yu.flatten()),uni.flatten(),\
                       (xf,yf))
    #Remove piston
    fig2 = fig2 - np.nanmean(fig2)
    uni2 = uni2 - np.nanmean(uni2)
    #Remove azimuthal variations to fairly high order
    fig2 = fig2 - fit.legendre2d(fig2,xo=4,yo=1)[0]
    uni2 = uni2 - fit.legendre2d(uni2,xo=4,yo=1)[0]
    #Determine coefficient needed to match axial sags
    fun = lambda c: rms(fig2-c*uni2)
    coeff = np.linspace(-5,5,1000)
    fom = [fun(c) for c in coeff]
    pdb.set_trace()
