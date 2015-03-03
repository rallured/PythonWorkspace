from numpy import *
from matplotlib.pyplot import *
import pdb
from conicsolve import primrad
from conicsolve import woltparam,primfocus
from plotting import myhist
import scattering as sc

#Compute dimple shape as a function of height and radius
def dimple(r,h=1.e-3,R=1.):
    ind = abs(r) > R
    r[ind] = R
    return (h/3)*(3-2*(r/R)**2-(r/R)**4+8*(r/R)**2*log(abs(r)/R))

def dimplep(r,h=1.e-3,R=1.):
    return (h/3)*(4*(r/R**2)-4*(r**3/R**4)+16*(r/R**2)*log(r/R))

def dimplepp(r,h=1.e-3,R=1.):
    return (h/3)*(20/R**2-12*r**2/R**4+16*log(r/R)/R**2)

#Compute radius of curvature for a dimple as a function of r
def dimplerad(r,h=1.e-3,R=1.):
    return abs((1+dimplep(r,h=h,R=R)**2)**(1.5)/dimplepp(r,h=h,R=R))

#Can we do this simple analysis solely in Python?
def dimpleWaveAnalysis(x=None,z=None,x0=None,h=14.e-3,R=12.):
    #Construct profile vectors
    if x is None:
        z = linspace(8400.,8500.,10000)
        x = primrad(z,220.,8400.) + dimple(z-8450.,h=h,R=R)

    #Compute other constants
    alpha,p,d,e = woltparam(220.,8400.)
    L1 = 100. #length of mirror
    wave = 1.e-6 #wavelength of light
    R0 = 220.
    DR = L1*sin(alpha) #aperture approximation
    dz = z[1]-z[0]

    #Construct psf vector
    psf = []
    if x0 is None:
        x0 = arange(-1.,1.,.001)
    dx0 = x0[1]-x0[0]
    for xf in x0:
        d2 = sqrt((x-xf)**2 + (z-8400.+primfocus(220.,8400.))**2)
        integrand = exp(-2j*pi*(d2-z)/wave)*dz*sqrt(x/d2)
        psf.append(abs(sum(integrand))**2 * DR/L1**2/wave/R0)
    psf = array(psf)

    #Compute standard deviation of PSF
    sigma = sum(abs(x0)*psf*dx0)
    
    return array(psf), sigma

#Treating dimple as PSD yield similar result?
def dimplePSDAnalysis(x=None,z=None,x0=None,h=14.e-3,R=12.):
    #Construct profile vectors
    if x is None:
        z = linspace(8400.,8500.,10000)
        dim = dimple(z-8450.,h=h,R=R)
        pwr = fft.fft(dim)/size(dim)
        pwr[0] = 0.
        pwr = abs(pwr)**2
        pwr = pwr[:floor(size(dim)/2+1)]*2
        freq = fft.fftfreq(size(dim),d=z[1]-z[0])
        freq = freq[:floor(size(dim)/2+1)]
        newdim = sc.invertPSD(z,freq,pwr)
        x = primrad(z,220.,8400.) + newdim

    #Compute other constants
    alpha,p,d,e = woltparam(220.,8400.)
    L1 = 100. #length of mirror
    wave = 1.e-6 #wavelength of light
    R0 = 220.
    DR = L1*sin(alpha) #aperture approximation
    dz = z[1]-z[0]

    #Construct psf vector
    psf = []
    if x0 is None:
        x0 = arange(-1.,1.,.001)
    dx0 = x0[1]-x0[0]
    for xf in x0:
        d2 = sqrt((x-xf)**2 + (z-8400.+primfocus(220.,8400.))**2)
        integrand = exp(-2j*pi*(d2-z)/wave)*dz*sqrt(x/d2)
        psf.append(abs(sum(integrand))**2 * DR/L1**2/wave/R0)
    psf = array(psf)

    #Compute standard deviation of PSF
    sigma = sum(abs(x0)*psf*dx0)
    
    return array(psf), sigma

#Compute slope error vector for perturbed vs nominal parabola
def slopeError(h=14.e-3,R=12.):
    z = linspace(8400.,8500.,1000000)
    nom = primrad(z,220.,8400.)
    pert = nom + dimple(z-8450.,h=h,R=R)
    return (diff(nom)-diff(pert))/(z[1]-z[0])

#Investigate fraction of light scattered beyond 1 arcsec
def scatterFractions(h=14.e-3,R=12.):
    #Compute wave optics scattered light
    psf,sig = dimpleWaveAnalysis(h=h,R=R)
    xo = arange(-1.,1.,.001)
    wavescat = sum(psf[abs(xo)<primfocus(220.,8400.)*5.e-6])*.001

    #Compute geometric convolution method
    psf,sig = dimpleWaveAnalysis(h=0.)
    sl = slopeError(h=h,R=R)
    off = 2*sl*primfocus(220.,8400.)
    bins = arange(round(min(off),3)-.0005,round(max(off),3)+.0015,.001)
    res = myhist(off,bins)
    dist = res[0].astype('float')/sum(res[0])
    geopsf = convolve(psf,dist)
    x = arange(round(min(off),3)-1,round(max(off),3)+1.000,.001)
    geoscat = sum(geopsf[abs(x)<primfocus(220.,8400.)*5.e-6])*.001

    return wavescat, geoscat
