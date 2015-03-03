from numpy import *
from matplotlib.pyplot import *
import pyfits,os,pdb
from zernikemod import stripnans
from utilities import fourier

#This module is for evaluating a multi-dimensional AXRO correction MTF
#Set a distortion map to axial sinusoids of varying amplitude, phase,
#and frequency

solvedir = '/Users/ryanallured/GoogleDrive/Optimization/optimization/solve_pzt/'
distortdir = '/Users/ryanallured/GoogleDrive/Optimization/example_data/'

#This function creates a distortion fits file and saves it to
#solve_pzt directory
def createDistortion(amp,freq,phase,filename):
    #Create distortion array
    y,x = mgrid[0:411,0:821]
    x,y = x*.5,y*.5
    d = amp*sin(2*pi*freq*y+phase)

    #Save as fits file in solve directory
    pyfits.writeto(distortdir+filename,d,clobber=True)

    return d

#Need to set distortion file, run solver, get resulting image,
#compute and return PSD
def correctedPSD(amp,freq,phase):
    #Set distortion array
    d = createDistortion(amp,freq,phase,'TestSine.fits')
    shade = pyfits.getdata(distortdir+'shade_mask.fits')
    d[shade==0] = nan
    d = stripnans(d)

    #Run solver to create solution
    os.chdir(solvedir)
    os.system(solvedir+'solve_pzt @fit-5x10mm-gap02-flex.par '
              'premath=Math.dat')

    #Load solution and ignore masked region
    resid = pyfits.getdata(solvedir+'X-resid.fits')
    resid[shade==0] = nan
    resid = stripnans(resid)

    #Compute average axial PSD
##    #Try removing piston from each slice first
##    for i in range(shape(resid)[0]):
##        resid[i] = resid[i] - mean(resid[i])
##    resid = resid - mean(resid) #Get rid of zero frequency component
    
##    fft2 = fft.fft2(resid*h2d)/size(resid) #Fourier components
##    axfft = 2*fft2[:,0] #Axial Fourier components, zero azimuthal frequency
##    axpsd = abs(axfft)**2
##    f = fft.fftfreq(size(axpsd),d=.5)
##    axpsd = axpsd[:size(axpsd)/2+1]
##    axpsd[0] = 0.
##    f = f[:size(f)/2+1]

    #Get windowed PSDs
    f,axpsdw = fourier.realPSD(resid,win=hanning,dx=.5)
    f = f[0] #Select only axial frequencies
    w = 2*pi*f/1000.
    axpsdw = axpsdw[:,0] #Select only axial frequencies

    f,origpsdw = fourier.realPSD(d,win=hanning,dx=.5)
    f = f[0]
    origpsdw = origpsdw[:,0]
    
    ripple = sum((w**2*axpsdw)[f>.15])/sum(w**2*origpsdw)

    #Get unwindowed PSDs for input reduction
    f,axpsd = fourier.realPSD(resid,dx=.5)
    f = f[0] #Select only axial frequencies
    axpsd = axpsd[:,0] #Select only axial frequencies

    f,origpsd = fourier.realPSD(d,dx=.5)
    f = f[0]
    origpsd = origpsd[:,0]
    
    correction = sum((w**2*axpsd)[f<.15])/sum(w**2*origpsd)
    
    return correction,ripple

def ampScan(amp,freq,phase):
    """Scan through frequencies and phases for a given amplitude
    amp is a scalar, freq and phase are arrays
    """
    corres = zeros((size(freq),size(phase)))
    ripres = copy(corres)
    for f in freq:
        for p in phase:
            corres[freq==f,phase==p],ripres[freq==f,phase==p] = \
                    correctedPSD(amp,f,p)
            sys.stdout.write('Freq: %.2e, Phase: %.2f' % (f,p))

    return corres,ripres
