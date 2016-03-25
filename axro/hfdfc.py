import numpy as np
import matplotlib.pyplot as plt
import pyfits
import axro.solver as slv
import pdb
import utilities.imaging.fitting as fitting
import utilities.imaging.man as man
import utilities.imaging.analysis as anal

mask = pyfits.getdata('/home/rallured/Dropbox/WFS/SystemAlignment/'
                      'HFDFC/160312_IFTest/Mask.fits')
ifs = pyfits.getdata('/home/rallured/Dropbox/WFS/SystemAlignment/'
                     'HFDFC/160312_IFTest/IFMeas/160313_MaskedIFs.fits')
ifs2 = pyfits.getdata('/home/rallured/Dropbox/WFS/SystemAlignment/'
                     'HFDFC/160312_IFTest/IFMeas/160313_MaskedIFs_NoTilt.fits')
ifs3 = pyfits.getdata('/home/rallured/Dropbox/WFS/SystemAlignment/HFDFC/160316IFs/160317_MedianSGSmooth.fits')

def maskandstrip(img):
    """Apply HFDFC mask and strip NaNs from an image"""
    img[mask==0] = np.nan
    return man.stripnans(img)

def createDist(amp,freq,phase,mag=9.75):
    """Creates a sinusoidal ripple distortion
    map for HFDFC.
    Frequency and amplitude need be in same units
    """
    #Create distortion array
    y,x = np.mgrid[0:128,0:128]
    x,y = x*(100./92),y*(100./92)
    d = amp*2*np.pi*freq*np.sin(2*np.pi*freq*y+phase)*mag

    #Apply mask
    d = maskandstrip(d)

    return d

def mapCTF(amp,freq,phase):
    """Loop through amplitude, frequency, and phase and
    compute CTF. This is brute force RMS slope before
    and after correction.
    """
    return [computeCTF(a,f,p) for a in amp for f in freq for p in phase]

def computeCTF(amp,freq,phase):
    """For a given amplitude, frequency, and phase of
    sinusoid, compute ratio of RMS slope before and after
    correction."""
    #Create distortion
    dist = createDist(amp,freq,phase)
    #Remove tilt
    dist = dist - np.mean(dist)
    #Compute pre-correction RMS slope
    prerms = anal.rms(dist*10000)

    #Apply correction
    cor = slv.rawOptimizer(ifs2*10000,dist*10000)

    #Compute post-correction RMS slope
    resid = dist*10000 - cor[0]
    postrms = anal.rms(resid)

    return postrms/prerms
