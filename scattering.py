from numpy import *
from matplotlib.pyplot import *
#from plotting import mycontour
import pdb
from conicsolve import primrad,primfocus,woltparam
import math

#Given an axial mirror profile and nominal Wolter Parameters,
#Compute PSF of focus from primary mirror only
#x is mirror radius at position z
#x0 is array of focal positions
def primaryPSF(x=None,z=None,perturb=None,x0=None,R0=220.,Z0=8400.,L1 =100.):
    #Construct profile vectors
    if x is None:
        z = linspace(8400.,8400.+L1,size(perturb))
        x = primrad(z,R0,Z0)
    if perturb is not None:
        x = x + perturb

    #Compute other constants
    alpha,p,d,e = woltparam(R0,Z0)
    wave = 1.e-6 #wavelength of light
    DR = L1*sin(alpha) #aperture approximation
    dz = z[1]-z[0]

    #Construct psf vector
    psf = []
    if x0 is None:
        x0 = arange(-1.,1.,.001)
    dx0 = x0[1]-x0[0]
    for xf in x0:
        d2 = sqrt((x-xf)**2 + (z-8400.+primfocus(R0,8400.))**2)
        integrand = exp(-2j*pi*(d2-z)/wave)*dz*sqrt(x/d2)
        psf.append(abs(sum(integrand))**2 * DR/L1**2/wave/R0)
    psf = array(psf)

    #Compute standard deviation of PSF
    sigma = sum(abs(x0)*psf*dx0)
    
    return array(psf)

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
