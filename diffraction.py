import numpy as np
import matplotlib.pyplot as plt
from numpy import pi,sqrt,exp
import pdb
from gaussfitter import twodgaussian

def quadPhase(wave,xg,yg,z):
    """Quadratic phase factor for meshgrid arrays xg,yg
    This is exp(1.j*k*r^2/2z)
    """
    return exp(1.j*(2*pi/wave)*(xg**2+yg**2)/(2*z))

def transferFn(z,wave,x,y):
    """Create transfer function for propagation
    convolution.
    This is -z*exp(ikR)/(2pi*R^2)*(ik-1/R)
    where R=sqrt((x-x0)^2+(y-y0)^2+z^2)
    Due to convolution integral, just take R(x0,y0,z)
    """
    R = sqrt(x**2+y**2+z**2)
    k = 2*pi/wave
    h = -z*exp(k*1j*R)/(2*pi*R**2)*(1j*k-1/R)
    return h

def rscProp(u0,wave,xg,yg,dz,dx=None):
    #Get spatial sampling
    if dx is None:
        dx = np.max(np.diff(xg.flatten()))
    #FFT of field distribution
    U0 = np.fft.fftn(u0)
    #FFT of transfer function
    tf = transferFn(dz,wave,xg,yg)
    H = np.fft.fftn(transferFn(dz,wave,xg,yg))
    #Inverse of convolution
    uz = np.fft.fftshift(np.fft.ifftn(U0*H))
    return uz

def propKernel(z,wave,fx,fy):
    """Create propagation kernel array for ASM propagation"""
    kern = exp(1.j*2*pi/wave*z*sqrt(1 - wave**2*(fx**2 + fy**2)))
    return kern

def asmProp(u0,wave,xg,yg,dz,dx=None):
    """Use angular spectrum propagation to propagate a
    field distribtion"""
    #FFT of field distribution
    U0 = np.fft.fftn(u0)
    #Get spatial sampling
    if dx is None:
        dx = np.max(np.diff(xg.flatten()))
    #Form frequency arrays
    sh = np.shape(xg)
    fx = np.fft.fftfreq(sh[0],d=dx)
    fy = np.fft.fftfreq(sh[1],d=dx)
    fy,fx = np.meshgrid(fy,fx)
    #Get propagation kernel
    kern = propKernel(dz,wave,fx,fy)
    #Inverse FFT for propagated field
    uz = np.fft.ifftn(U0*kern)
    return uz

def cdaMask(gridwidth,N,subx,suby):
    """Determine amount of power making it through the
    CDA double pass system.
    """
    #Define Gaussian function for beam
    fn = twodgaussian([0.,1.,0.,0.,.75,.75,0.])
    #Define pencil beam 1/e^2 diameter 3 mm
    g = np.linspace(-gridwidth,gridwidth,N)
    x,y = np.meshgrid(g,g)
    beam = fn(x,y)
    beam = beam / np.sum(beam)
    beam = sqrt(beam)
    #Subaperture .6 mm circular pupil
    beam[np.logical_or(np.abs(x)>subx,\
                       np.abs(y)>suby)] = 0.
    pdb.set_trace()
    #Propagate 370 mm * 2 back to aperture
    beam = asmProp(beam,632.e-6,x,y,370.*2)
    #Subaperture .6 mm circular pupil
    beam[np.logical_or(np.abs(x)>subx,\
                       np.abs(y)>suby)] = 0.
    #Return power transmitted through mask upon second pass
    return np.sum(np.abs(beam)**2)

##Test out simulating lens focusing a perfect plane wave
##wave = 633.e-6 #633 nm
##k = 2*pi/wave
##
##Define Gaussian beam, make sure phase array is larger than beam
##This prevents aliasing
##meshhalfwidth = 10.
##N = 2000
##x = np.linspace(-meshhalfwidth,meshhalfwidth,N)
##y = np.copy(x)
##
##xg,yg = np.meshgrid(x,y)
##beamwidth = 2.
##beam = exp(-(xg**2+yg**2)/(2*pi*beamwidth))
##
##Cylindrical beam
##ind = np.where(sqrt(xg**2+yg**2)<beamwidth)
##beam = np.zeros(np.shape(xg))
##beam[ind] = 1.
##
##Multiply by quadratic phase factor for focusing lens
##f = 500.
##beam = beam * quadPhase(wave,xg,yg,-500.)
##
##Propagate with RSC
##newph = rscProp(beam,wave,xg,yg,500.)
##newph2 = asmProp(beam,wave,xg,yg,500.)

#Propagate to focus, multiply by quadratic phase factor
#and take FFT
#What does frequency mean? There will be scaling to do
#We're at the focus, so we know roughly what the spot size should be
##d = 500.
##beam = beam * exp(1.j*k/(2*d)*(xg**2+yg**2)) #Basically cancels out the lens
##beamfoc = np.fft.fft2(beam)
##freq = np.fft.fftshift(np.fft.fftfreq(np.size(x),x[1]-x[0])*wave*d)
##
##beammag = np.fft.fftshift(np.real(beamfoc*np.conjugate(beamfoc)))
