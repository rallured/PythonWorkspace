import numpy as np
import matplotlib.pyplot as plt
#from plotting import mycontour
import pdb
from traces.conicsolve import primrad,primfocus,woltparam
import math
import utilities.imaging.man as man
import scatter
import inducedPolarization as ip
import time
import utilities.fourier as fourier

import numba
#from numba import cuda

def computeHEW(x0,psf):
    """
    Compute the HEW from a PSF via the CDF method
    """
    dx2 = np.diff(x0)[0]
    cdf = np.cumsum(psf)*dx2
    hpdPSF = x0[np.argmin(np.abs(cdf-.75))]-\
             x0[np.argmin(np.abs(cdf-.25))]
    return hpdPSF

def primary2DPSF(img,dx,R0=220.,Z0=8400.,x0=np.linspace(-5.,5.,1001),\
                 wave=1.24e-6):
    """
    Create height vector and radius img based on
    radial distortion input data.
    Then pass to a F2PY scattering function to
    compute PSF over observation points.
    """
    tstart = time.time()
    #Remove NaNs if they exist
    img = man.stripnans(img)
    #Create height vector
    graze = woltparam(R0,Z0)[0]
    foc = primfocus(R0,Z0)
    z = np.arange(np.shape(img)[0])*dx*np.cos(graze)+Z0
    #Create radial position img
    rad = primrad(z,R0,Z0)
    rad2 = np.flipud(np.transpose(np.tile(rad,(np.shape(img)[1],1))))
    distortion = np.transpose(rad2-img/1e3)
    z = z[::-1]
    #Compute length for each slice
    length = np.array([(np.sum(~np.isnan(li))-1)*dx \
                       for li in distortion],order='F')
    DR = length*np.sin(graze)
    #Integrate each slice in Fortran
    pdb.set_trace()
    psf = scatter.primarypsf(distortion,z-Z0,length,x0,wave,foc,R0,graze)
    print time.time()-tstart
    
    return psf
    
##def primaryGPU2D(img,dx,R0=220.,Z0=8400.,x0=np.linspace(-5.,5.,1001)):
##    """
##    Create height vector and radius img based on
##    radial distortion input data.
##    Then pass to a F2PY scattering function to
##    compute PSF over observation points.
##    """
##    tstart = time.time()
##    #Remove NaNs if they exist
##    img = man.stripnans(img)
##    #Create height vector
##    graze = woltparam(R0,Z0)[0]
##    foc = primfocus(R0,Z0)
##    z = np.arange(np.shape(img)[0])*dx*np.cos(graze)+Z0
##    #Create radial position img
##    rad = primrad(z,R0,Z0)
##    rad2 = np.flipud(np.transpose(np.tile(rad,(np.shape(img)[1],1))))
##    distortion = np.transpose(rad2-img/1e3)
##    z = z[::-1]
##    #Compute length for each slice
##    length = np.array([(np.sum(~np.isnan(li))-1)*dx \
##                       for li in distortion],order='F')
##    DR = length*np.sin(graze)
##    psf = [gpuPSF(np.ascontiguousarray(z),\
##                  np.ascontiguousarray(ri),x0,wave=1.24e-6,Z0=8400.,R0=220.) \
##           for ri in distortion]
##
##    print time.time()-tstart
##    
##    return psf
    

#Given an axial mirror profile and nominal Wolter Parameters,
#Compute PSF of focus from primary mirror only
#x is mirror radius at position z
#x0 is array of focal positions
def primaryPSF(x=None,z=None,perturb=None,x0=np.arange(-1.,1.,.001)\
               ,R0=220.,Z0=8400.,wave=1.24e-6):
    #Construct profile vectors
    if z is None:
        z = np.arange(np.shape(img)[0])*dx*np.cos(graze)+Z0
    if x is None:
        x = primrad(z,R0,Z0)
    if perturb is not None:
        pdb.set_trace()
        x = x + perturb

    L1 = z.max()-z.min()

    #Get rid of NaNs
    ind = ~np.isnan(x)
    #Rescale length based on amount of missing data
    L1 = L1 * float(np.sum(ind))/len(z)
    x = x[ind]
    z = z[ind]

    #Compute other constants
    alpha,p,d,e = woltparam(R0,Z0)
    DR = L1*np.sin(alpha) #aperture approximation
    dz = z[1]-z[0]

    #Construct psf vector
    psf = []
    if x0 is None:
        x0 = np.arange(-1.,1.,.001)
    dx0 = x0[1]-x0[0]
    pdb.set_trace()
    for xf in x0:
        d2 = np.sqrt((x-xf)**2 + (z-Z0+primfocus(R0,Z0))**2)
        integrand = np.exp(-2j*np.pi*(d2-(z-Z0))/wave)*dz*np.sqrt(x/d2)
        psf.append(abs(np.sum(integrand))**2 * DR/L1**2/wave/R0)
    psf = np.array(psf)

    #Compute standard deviation of PSF
    sigma = np.sum(abs(x0)*psf*dx0)
    
    return np.array(psf)

##@cuda.jit('void(double[:],double[:],double[:],double[:],double)')
##def integratePSF(zin,rin,x0,outPSF,wave):
##    """
##    Divide a PSF computation amongst GPU threads by observation point.
##    Each thread does its own reduction, eliminating interthread communication.
##    Assume data has been formulated properly, and is in zin,rin,x0 format.
##    Only computation to be done is Raimondi integral.
##    Multiplicative constants handled outside the integral in wrapper function.
##    Shared memory technique did not work due to limited block memory.
##    Unclear whether restructuring computation to load distortion
##    in shared memory will help.
##    """
####    #Declare shared memory
####    zarr = cuda.shared.array(shape=100000,dtype=numba.double)
####    rarr = cuda.shared.array(shape=100000,dtype=numba.double)
####    #Get distortion data into shared memory
####    i = cuda.grid(1)
####    if i < len(zin):
####        zarr[i] = zin[i]
####        rarr[i] = rin[i]
####    else:
####        zarr[i] = 0.
####        rarr[i] = 0.
####    cuda.syncthreads()
##    
##    #Grab observation point for this thread
##    i = cuda.grid(1)
##    thisx = x0[i]
##
##    #Loop through distortion data and compute PSF
##    integrand = np.complex128(0.)
##    for n in range(len(zin)):
##        if rin[n] > 0.:
##            d2 = math.sqrt((thisx-rin[n])**2+zin[n]**2)
##            expon = (-2*np.pi/wave)*(d2-zin[n])
##            integrand += math.sqrt(rin[n]/d2) * \
##                        (math.cos(expon)+1j*math.sin(expon))
##            #outPSF[i] = d2
##            
##
##    #Compute PSF
##    if i < len(x0):
##        outPSF[i] = abs(integrand)**2#*dz**2*dR/L**2/wave/R0
##    
##    return None
##
##def gpuPSF(z,r,x0,wave=1.24e-6,Z0=8400.,R0=220.):
##    """
##    Wrapper for GPU PSF calculation
##    """
####    #Testing code
####    z = np.linspace(8400.,8500.,10000)
####    r = primrad(z,220.,8400.)
####    z = z-8400.+primfocus(220.,8400.)
####    x0 = np.linspace(-1.,1.,10000)
##
##    #Get contant parameters
##    print 'Tstart: %.4f' % time.time()
##    dz = z[1]-z[0]
##    alpha = .25*np.arctan(R0/Z0)
##    L = np.sum(~np.isnan(r))*dz
##    dR = L*np.sin(alpha)
##    
##    psf = np.zeros(len(x0))
##    psfg = cuda.to_device(psf)
##    zg = cuda.to_device(z)
##    rg = cuda.to_device(r)
##    x0g = cuda.to_device(x0)
##
##    stream = cuda.stream()
##    tstart = time.time()
##    integratePSF[len(x0)/64+1,64](zg,rg,x0g,psfg,wave)
##
##    psfg.copy_to_host(psf)
##
##    print 'Tend: %.4f' % time.time()
##
##    return psf*dz**2*dR/L**2/wave/R0


#Return AXRO PSD requirement
#f is frequency vector in 1/mm units
def axroPSD(f):
    sigma_1 = 1.5 * 10**-4
    L1 = 35

    sigma_2 = 4.4 * 10**-6
    L2 = 20

    PSD = 2 * sigma_1**2  * L1 * exp(-pi * f * L1)**2 +\
          4 * sigma_2**2 * L2 / (1 + (2 * pi * f * L2)**2 )

    return PSD

#Return integral of AXRO PSD from f1 to f2 (in 1/mm)
def intAxroPsd(f1,f2):
    sigma_1 = 1.5 * 10**-4
    L1 = 35.

    sigma_2 = 4.4 * 10**-6
    L2 = 20.

    gausterm = sigma_1**2*L1/sqrt(pi) * \
               (math.erf(pi*L1*f2)-math.erf(pi*L1*f1))
    lorentzterm = 2*sigma_2**2/pi * (arctan(2*pi*L2*f2)-arctan(2*pi*L2*f1))

    return gausterm+lorentzterm

#Invert PSD (in units of mm^2) to a random surface
#Requires frequency vector in units of 1/mm
#Adds sinusoid of frequency f_i, magnitude sqrt(psd_i),
#and a random phase between 0 and 2*pi
def invertPSD(x,freq,psd):
    #Add in cosines of appropriate magnitudes, frequencies, with random phase
    surf = zeros(size(x))
    for i in range(size(freq)):
        surf = surf + sqrt(psd[i])*cos(2*pi*freq[i]*x+random.rand()*2*pi)
    return surf

#Calculate scattering coefficients given a 1D sinusoidal surface
#Takes incidence angle (t1 in degrees),
#wavelength of light, period of sinusoid, and amplitude (h) of sinusoid
#as parameters
#Outputs diffracted orders, scattering coefficients, and angular half-widths
#This solution is given in sec. 4.3 of Beckmann and Spizzichino (1963)
#Also in sec. 2.7 of Maystre's review of grating theory in Progress
#in Optics (1984)
def sinscatter(t1,wave,period,h):
    #Calculate intermediate parameters
    t1 = t1 * pi/180
    m = arange(-100,100)
    t2m = m*wave/period+sin(t1) #sine of diffraction angle for order m
    
    m = m[where(abs(t2m)<1)] #remove evanescent orders
    t2m = t2m[where(abs(t2m)<1)] #remove evanescent orders
    t2m = arcsin(t2m)
    
    s = (2*pi/wave)*h*(cos(t1)+cos(t2m))
    rho = (1/cos(t2m))*(1+cos(t1+t2m))/(cos(t1)+cos(t2m))*((-1j)**m)*sp.jv(m,s)

    print wave**2*h/period**3

    

    return m,rho*conjugate(rho)*cos(t2m)/cos(t1),t2m*180/pi

#Input initial angle of incidence (t1), rms roughness amplitude (sigma),
#correlation length (T), wavelength (wave)
#Use consistent units for all spatial parameters
#Outputs mean scattered power as a function of scattering angle
#View in a contour plot, should be sharply peaked about specular
#t2 and t3 are maximum angular deviations to calculate over
#(i.e. from -t2 to +t2)
#You can change the 1000 to something else for more coarse or fine
#calculation grid
def normalscatter(wave,t1,sigma,T,t2,t3):
    #Initialize vectors of scattering angles
    t1 = t1 * pi/180.
    t2 = t2 * pi/180.
    t3 = t3 * pi/180.
    t2 = linspace(t1-t2,t1+t2,1000)
    t3 = linspace(-t3,t3,1000)

    #Create 2D arrays where second index is in plane scatter angle
    t2,t3 = meshgrid(t2,t3)

    F = lambda t1,t2,t3: (1.+cos(t1)*cos(t2)-sin(t1)*sin(t2)*cos(t3)) /\
        (cos(t1)*(cos(t1)+cos(t2)))

    #Compute factors in scattered power expression
    lx = 1.
    ly = 1.
    A = lx*ly #1 square meter
    k = 2*pi/wave
##    F = (1.+cos(t1)*cos(t2)-sin(t1)*sin(t2)*cos(t3)) /\
##        (cos(t1)*(cos(t1)+cos(t2)))
    vz = -k*(cos(t1)+cos(t2))
    vy = -k*(sin(t2)*sin(t3))
    vx = k*(sin(t1)-sin(t2)*cos(t3))
    vxy = sqrt(vx**2+vy**2)
##    vxy = k*sqrt(sin(t1)**2-2*sin(t1)*sin(t2)*cos(t3)+sin(t2)**2)
    g = sqrt((vz*sigma)**2)

    #Put it together into scattered power
    #Default to rough surface formula
    fact = pi*F(t1,t2,t3)**2*T**2 / (A*vz**2*sigma**2)
    expterm = exp(-vxy**2*T**2/(4*(vz**2)*sigma**2))
    power = fact*expterm
    #If g << 1, it makes sense to use smooth surface formula
    ind = where(g < 1)
    rho0 = sinc(vx*lx)*sinc(vy*ly)
    smooth = exp(-g) * (rho0**2 + pi*T**2*F(t1,t2,t3)/A * exp(-vxy**2*T**2/4.))
    power[ind] = smooth[ind]

    #Comment this out if you don't want to plot the distribution
    clf()
    contourf(power,levels=linspace(0.,power.max(),100))
    colorbar(format='%.3e')

    return power,t2,t3,fact,expterm

#Using simplified theory by Aschenbach for separating geometric
#vs. physical optics
def septest(t1,t2):
    #Define power spectrum
    #frequency in inverse microns
    #roughness power in nm^3
    psd = lambda freq: .5/(freq**2.2)

    #Compute discrete power spectrum
    t2 = linspace(t1-t2,t1+t2,10000) * pi/180
    t1 = t1 * pi/180
    freq = 2*pi/1.e-3*(sin(t2)-sin(t1))
    power = psd(abs(freq))
    #Set all power above 5 mm period to 0.
    ind = where(abs(freq) < 1./5.e3)
    power[ind] = 0.
    powsum = 0.
    for i in arange(1,size(power)):
        powsum += power[i]*abs(freq[i]-freq[i-1])/1000.
    print sqrt(powsum)
    
    #Compute scatter intensity vs scatter angle theta 2
    scatter = 4*(2*pi/1.)**3*cos(t1)*cos(t2)**2*power

    pdb.set_trace()

    return t2,scatter,power

def computeAperture(R0=220.,L=100.,Z0=1e4,wave=1e-6):
    x0 = np.arange(-.2,.2,.0001)
    return computeHEW(x0,primaryPSF(z=np.linspace(Z0,Z0+L,100),Z0=Z0,\
                             R0=R0,wave=wave,x0=x0))/primfocus(R0,Z0)\
                             *180/np.pi*60**2

def investigateAperture(L=100.):
    """
    Run through XRS energies and radii and compute aperture diffraction.
    """
    x0 = np.arange(-.2,.2,.0001)
    energy = np.linspace(.3,10.,10)

    wave = 1240.e-9/energy
    radii = np.linspace(300.,1000.,10)
    z = np.linspace(1e4,1e4+L,100)
    
    hew = [[computeHEW(x0,primaryPSF(z=z,Z0=1e4,R0=r,wave=w,x0=x0))/\
          primfocus(r,1e4)*180/np.pi*60**2 for r in radii] for w in wave]

    plt.figure()
    [plt.semilogy(radii,h,label='%.2f' % e) for h,e in zip(hew,energy)]
    plt.legend(loc='lower left')

    return hew

def weightedAperture(energy=1000.,L=100.):
    """

    """
    #Get Rx
    rx = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO/WSTracing/'
                                    '150528_Pauls_Rx.csv',delimiter=','))
    geo = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO/WSTracing/'
                                     'geometric_transmission_102711.txt'))
    rx = rx[:,geo[1]>0]
    geo = geo[1][geo[1]>0]

    #Compute weights
    nk = ip.readNK('/home/rallured/Dropbox/inducedPol/Ir_llnl_cxro.nk')
    ang = .25*np.arctan(rx[1]/1e4)
    refl = ip.computeReflectivities(ang,energy,.5,1.,nk)[0]
    weights = rx[-1]*refl**2

    #Compute shell diffraction
    app = np.array([computeAperture(R0=r,L=L,Z0=1e4,wave=1240.e-6/energy)\
                    for r in rx[1]])
    Z0 = 1e4

    x0 = np.arange(-.2,.2,.0001)
    psfs = np.array([primaryPSF(z=np.linspace(Z0,Z0+L,100),Z0=Z0,\
                             R0=r,wave=1240.e-6/energy,x0=x0)/primfocus(r,Z0)\
                             *180/np.pi*60**2 for r in rx[1]])
    psf = np.average(psfs,axis=0,weights=weights)

    #Return weighted sum
    return np.average(app,weights=weights),psf

def roughnessHPD(freq,psd,R0=220.,Z0=10000.,x0=np.linspace(-5.,5.,1001),\
                 wave=1.24e-6):
    """
    Compute the HPD of a random profile due to roughness
    """
    #Construct random profile
    prof = np.zeros((len(freq),5))
    for i in range(5):
        prof[:,i] = fourier.randomProfile(freq,psd)

##    prof = fourier.randomProfile(freq,psd)
##    prof = prof.reshape((len(freq),1))

    #Get PSF
    dx = 1./max(abs(freq))/2.
    print dx
    psf = primary2DPSF(prof*1e3,dx,R0=R0,Z0=Z0,x0=x0,wave=wave)

    return psf,prof
