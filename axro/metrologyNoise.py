import numpy as np
import matplotlib.pyplot as plt
import utilities.metrology as met
import utilities.fourier as fourier
import pdb,glob
import astropy.io.fits as pyfits
from utilities.imaging.fitting import legendre2d
import os

def flatNoiseCGH():
    """
    Quantify and visualize impact of parroting for
    Solar-B flat measurements (through CGH)
    """
    #Get data
    wdir = '/home/rallured/Dropbox/AXRO/Metrology/NoiseStudy/FlatMeasurements/'
    d1,dx1 = met.read4DFits(wdir+'161205_RefFlat_Avg8_Meas1.fits')
    d2,dx2 = met.read4DFits(wdir+'161205_RefFlat_Avg8_Meas2.fits')
    p1,px1 = met.read4DFits(wdir+'161205_RefFlat_ParrotingTestPitch_Meas1.fits')
    p2,px2 = met.read4DFits(wdir+'161205_RefFlat_ParrotingTestPitch_Meas2.fits')
    p3,px3 = met.read4DFits(wdir+'161205_RefFlat_ParrotingTestRoll_Meas1.fits')
    p4,px4 = met.read4DFits(wdir+'161205_RefFlat_ParrotingTestRoll_Meas2.fits')

    #Construct baseline power spectra
    f1,pow1 = fourier.meanPSD(d1-d2,win=np.hanning,dx=dx1)
    f2,pow2 = fourier.meanPSD(d1-d2,win=np.hanning,dx=dx1,axis=1)
    
    #Construct parroted power spectra
    f3,pow3 = fourier.meanPSD(p1-p2,win=np.hanning,dx=dx1)
    f4,pow4 = fourier.meanPSD(p1-p2,win=np.hanning,dx=dx2,axis=1)
    f5,pow5 = fourier.meanPSD(p3-p4,win=np.hanning,dx=dx1)
    f6,pow6 = fourier.meanPSD(p3-p4,win=np.hanning,dx=dx2,axis=1)

    #Plot
    plt.loglog(f1,pow1/f1[0],label='Axial Baseline')
    plt.loglog(f2,pow2/f2[0],label='Azimuthal Baseline')
    plt.loglog(f3,pow3/f3[0],label='Pitch Axial')
    plt.loglog(f4,pow4/f4[0],label='Pitch Azimuthal')
    plt.loglog(f5,pow5/f5[0],label='Roll Axial')
    plt.loglog(f6,pow6/f6[0],label='Roll Azimuthal')
    plt.title('Residual Fringe Repeatability Impact')
    plt.xlabel('Frequency (1/mm)')
    plt.ylabel('Power ($\mu$m$^2$ mm)')
    plt.grid()
    plt.legend(loc='lower left')

    return f1,pow1

def refCylNoise():
    """
    Quantify and visualize impact of parroting for
    PCO cylinder measurements (through CGH)
    """
    #Get data
    wdir = '/home/rallured/Dropbox/AXRO/Metrology/' \
           'NoiseStudy/RefCylinderMeasurements/'
    d1,dx1 = met.read4DFits(wdir+'161205_RefCylinder_Avg8_Meas1.fits')
    d2,dx2 = met.read4DFits(wdir+'161205_RefCylinder_Avg8_Meas2.fits')
    d3,dx3 = met.read4DFits(wdir+'161205_RefCylinder_Avg8_Meas3.fits')

    p1,px1 = met.read4DFits(wdir+'161205_RefCylinder_'
                            'ParrotingTestPitch_Meas1.fits')
    p2,px2 = met.read4DFits(wdir+'161205_RefCylinder_'
                            'ParrotingTestPitch_Meas2.fits')
    p3,px3 = met.read4DFits(wdir+'161205_RefCylinder_'
                            'ParrotingTestRoll_Meas1.fits')
    p4,px4 = met.read4DFits(wdir+'161205_RefCylinder_'
                            'ParrotingTestRoll_Meas2.fits')

    #Construct baseline power spectra
    f1,pow1 = fourier.meanPSD(d1-d2,win=np.hanning,dx=dx1)
    f2,pow2 = fourier.meanPSD(d1-d2,win=np.hanning,dx=dx1,axis=1)
    
    #Construct parroted power spectra
    f3,pow3 = fourier.meanPSD(p1-p2,win=np.hanning,dx=dx1)
    f4,pow4 = fourier.meanPSD(p1-p2,win=np.hanning,dx=dx2,axis=1)
    f5,pow5 = fourier.meanPSD(p3-p4,win=np.hanning,dx=dx1)
    f6,pow6 = fourier.meanPSD(p3-p4,win=np.hanning,dx=dx2,axis=1)

    #Plot
    plt.loglog(f1,pow1/f1[0],label='Axial Baseline')
    plt.loglog(f2,pow2/f2[0],label='Azimuthal Baseline')
    plt.loglog(f3,pow3/f3[0],label='Pitch Axial')
    plt.loglog(f4,pow4/f4[0],label='Pitch Azimuthal')
    plt.loglog(f5,pow5/f5[0],label='Roll Axial')
    plt.loglog(f6,pow6/f6[0],label='Roll Azimuthal')

    return f1,pow1

def flatNoisePellicle():
    """
    Quantify and visualize noise through pellicle observing solar B flat
    """
    #Get data
    wdir = '/home/rallured/Dropbox/AXRO/Metrology/' \
           'NoiseStudy/SolarBwPellicle/'
    d1,dx1 = met.read4DFits(wdir+'161209_Avg8_Meas1.fits')
    d2,dx2 = met.read4DFits(wdir+'161209_Avg8_Meas2.fits')
    d3,dx3 = met.read4DFits(wdir+'161209_Avg8_Meas3.fits')
    d4,dx4 = met.read4DFits(wdir+'161209_Avg8_Meas4.fits')

    #Construct power spectra
    f12,pow12 = fourier.meanPSD((d1-d2)[:,100:-100],\
                                win=np.hanning,dx=dx1,irregular=True)
    f23,pow23 = fourier.meanPSD((d2-d3)[:,100:-100],\
                                win=np.hanning,dx=dx1,irregular=True)
    f34,pow34 = fourier.meanPSD((d3-d4)[:,100:-100],\
                                win=np.hanning,dx=dx1,irregular=True)
    f14,pow14 = fourier.meanPSD((d1-d4)[:,100:-100],\
                                win=np.hanning,dx=dx1,irregular=True)

    #Mid frequency
    midfreq = [1000*np.sqrt(np.sum(p[np.logical_and(f>.1,f<1.)])) \
               for f,p in zip([f12,f23,f34,f14],[pow12,pow23,pow34,pow14])]

    #Plot
    plt.loglog(f12,pow12/f12[0],label='1-2: %.2f' % midfreq[0])
    plt.loglog(f23,pow23/f23[0],label='2-3: %.2f' % midfreq[1])
    plt.loglog(f34,pow34/f34[0],label='3-4: %.2f' % midfreq[2])
    plt.loglog(f14,pow14/f14[0],label='1-4: %.2f' % midfreq[3])
    plt.legend(loc='lower left')
    plt.grid()
    plt.title('4D Repeatability: SolarB Flat+Pellicle')
    plt.xlabel('Frequency (1/mm)')
    plt.ylabel('Power ($\mu$m$^2$ mm)')

    print midfreq

    return f12,pow12

def PCO1S12Noise():
    """
    Quantify and visualize noise through pellicle observing solar B flat
    """
    #Get data
    wdir = '/home/rallured/Dropbox/AXRO/Metrology/' \
           'NoiseStudy/TestOptics_PCO1S12/'
    d1,dx1 = met.read4DFits(wdir+'161202_PCO1S12_4InchCut_Avg8_Meas1.fits')
    d2,dx2 = met.read4DFits(wdir+'161202_PCO1S12_4InchCut_Avg8_Meas2.fits')
    d3,dx3 = met.read4DFits(wdir+'161202_PCO1S12_4InchCut_Avg8_Meas3.fits')

    #Construct power spectra
    f12,pow12 = fourier.meanPSD((d1-d2)[:,100:-100],\
                                win=np.hanning,dx=dx1,irregular=True)
    f23,pow23 = fourier.meanPSD((d2-d3)[:,100:-100],\
                                win=np.hanning,dx=dx1,irregular=True)
    f13,pow13 = fourier.meanPSD((d1-d3)[:,100:-100],\
                                win=np.hanning,dx=dx1,irregular=True)

    #Mid frequency
    midfreq = [1000*np.sqrt(np.sum(p[np.logical_and(f>.1,f<1.)])) \
               for f,p in zip([f12,f23,f13],[pow12,pow23,pow13])]

    #Plot
    plt.loglog(f12,pow12/f12[0],label='1-2: %.2f' % midfreq[0])
    plt.loglog(f23,pow23/f23[0],label='2-3: %.2f' % midfreq[1])
    plt.loglog(f13,pow13/f13[0],label='1-3: %.2f' % midfreq[2])
    plt.legend(loc='lower left')
    plt.grid()
    plt.title('4D Repeatability: PCO1S12')
    plt.xlabel('Frequency (1/mm)')
    plt.ylabel('Power ($\mu$m$^2$ mm)')

    print midfreq

    return f12,pow12

def compareRepeatability():
    """
    Plot repeatability spectra of various configurations
    """
    fc,pc = refCylNoise()
    fp,pp = flatNoisePellicle()
    fcgh,pcgh = flatNoiseCGH()
    fpco,ppco = PCO1S12Noise()

    midfreq = [1000*np.sqrt(np.sum(p[np.logical_and(f>.1,f<1.)])) \
               for f,p in zip([fp,fc,fcgh,fpco],[pp,pc,pcgh,ppco])]
    plt.clf()
    
    plt.loglog(fp,pp/fp[0],label='Flat+Pellicle: %.2f' % midfreq[0])
    plt.loglog(fc,pc/fc[0],label='Cylinder+CGH: % .2f' % midfreq[1])
    plt.loglog(fcgh,pcgh/fcgh[0],label='Flat+CGH: %.2f' % midfreq[2])
    plt.loglog(fpco,ppco/fpco[0],label='PCO1S12: %.2f' % midfreq[3])
    plt.title('Repeatability Comparison')
    plt.xlabel('Frequency (1/mm)')
    plt.ylabel('Power ($\mu$m$^2$ mm)')
    plt.grid()
    plt.legend(loc='lower left')

def investigate4DRepeatability():
    """
    Analyze repeatability taken on the SolarB flat
    Temporal mode and Dynamic mode with fringes nulled and
    with an axial tilt.
    Plot figure PSDs with progressive averaging for each mode
    Plot comparative repeatability PSDs
    Plot comparative figure PSDs with 32 averages
    """
    parentdir = '/home/rallured/Dropbox/Interferometer/SolarBFlat/Repeatability/'
    avgs = [1,2,4,8,16,32]

    #Temporal with fringes tilted
    fn = glob.glob(parentdir+'Tilt/17*RepeatabilityTiltTemporal*.bin')
    fn.sort()
    dx = met.readFlatScript(fn[0].split('.')[0])[1]
    d = np.array([met.readFlatScript(fi.split('.')[0])[0] for fi in fn])
    #Make progressive averaging plot
    plt.figure('TemporalTiltedFigure')
    for i in np.arange(6)*2:
        f,p = fourier.meanPSD(d[i],win=np.hanning,dx=dx,irregular=True,\
                              minpx=200)
        plt.loglog(f,p/f[0],label=str(avgs[i/2]))
    plt.legend(loc='lower left')
    plt.title('Solar B PSD - Temporal,Tilted')
    plt.xlabel('Frequency (1/mm)')
    plt.ylabel('Power ($\mu$m$^2$ mm)')
    plt.grid()
    #Get repeatability
    reptemptilt = d[-1]-d[-2]
    figtemptilt = d[-1]

    #Dynamic with fringes tilted
    fn = glob.glob(parentdir+'Tilt/17*RepeatabilityTilt_*.bin')
    fn.sort()
    dx = met.readFlatScript(fn[0].split('.')[0])[1]
    d = [met.readFlatScript(fi.split('.')[0])[0] for fi in fn]
    #Make progressive averaging plot
    plt.figure('DynamicTiltedFigure')
    for i in np.arange(6)*2:
        f,p = fourier.meanPSD(d[i],win=np.hanning,dx=dx,irregular=True,\
                              minpx=200)
        plt.loglog(f,p/f[0],label=str(avgs[i/2]))
    plt.legend(loc='lower left')
    plt.title('Solar B PSD - Dynamic,Tilted')
    plt.xlabel('Frequency (1/mm)')
    plt.ylabel('Power ($\mu$m$^2$ mm)')
    plt.grid()
    #Get repeatability
    repdyntilt = d[-1]-d[-2]
    figdyntilt = d[-1]
    
    #Temporal with fringes nulled
    fn = glob.glob(parentdir+'Nulled/17*.bin')
    fn.sort()
    dx = met.readFlatScript(fn[0].split('.')[0])[1]
    d = np.array([met.readFlatScript(fi.split('.')[0])[0] for fi in fn])
    #Make progressive averaging plot
    plt.figure('TemporalNulledFigure')
    for i in np.arange(6)*2:
        f,p = fourier.meanPSD(d[i],win=np.hanning,dx=dx,irregular=True,\
                              minpx=200)
        plt.loglog(f,p/f[0],label=str(avgs[i/2]))
    plt.legend(loc='lower left')
    plt.title('Solar B PSD - Temporal,Nulled')
    plt.xlabel('Frequency (1/mm)')
    plt.ylabel('Power ($\mu$m$^2$ mm)')
    plt.grid()
    #Get repeatability
    reptempnull = d[-1]-d[-2]
    figtempnull = d[-1]
    
    #Dynamic with fringes nulled
    d = pyfits.getdata('/home/rallured/Dropbox/Interferometer/'
                       'SolarBFlat/Repeatability/'
                       'Nulled/170103_Processed.fits')
    rep = np.array([d[i,0]-d[i,1] for i in range(32)])
    #Make progressive averaging plot
    plt.figure('DynamicNulledFigure')
    for i in [0,1,3,7,15,31]:
        f,p = fourier.meanPSD(d[i,0],win=np.hanning,dx=dx,irregular=True,\
                              minpx=200)
        plt.loglog(f,p/f[0],label=str(i+1))
    plt.legend(loc='lower left')
    plt.title('Solar B PSD - Dynamic,Nulled')
    plt.xlabel('Frequency (1/mm)')
    plt.ylabel('Power ($\mu$m$^2$ mm)')
    plt.grid()
    #Get repeatability
    repdynnull = d[-1][0]-d[-1][1]
    figdynnull = d[-1][0]

    #Make comparative repeatability plots with 32 averages
    plt.figure('CompareRepeatability')
    f,p = fourier.meanPSD(repdynnull,win=np.hanning,dx=dx,irregular=True,\
                          minpx=200)
    plt.loglog(f,p/f[0],label='Dynamic,Nulled')
    f,p = fourier.meanPSD(repdyntilt,win=np.hanning,dx=dx,irregular=True,\
                          minpx=200)
    plt.loglog(f,p/f[0],label='Dynamic,Tilted')
    f,p = fourier.meanPSD(reptemptilt,win=np.hanning,dx=dx,irregular=True,\
                          minpx=200)
    plt.loglog(f,p/f[0],label='Temporal,Tilted')
    f,p = fourier.meanPSD(reptempnull,win=np.hanning,dx=dx,irregular=True,\
                          minpx=200)
    plt.loglog(f,p/f[0],label='Temporal,Nulled')
    plt.legend(loc='lower left')
    plt.title('Solar B Repeatability - 32 Averages')
    plt.xlabel('Frequency (1/mm)')
    plt.ylabel('Power ($\mu$m$^2$ mm)')
    plt.grid()

    #Make comparative figure plots with 32 averages
    plt.figure('CompareFigure')
    f,p = fourier.meanPSD(figdynnull,win=np.hanning,dx=dx,irregular=True,\
                          minpx=200)
    plt.loglog(f,p/f[0],label='Dynamic,Nulled')
    f,p = fourier.meanPSD(figdyntilt,win=np.hanning,dx=dx,irregular=True,\
                          minpx=200)
    plt.loglog(f,p/f[0],label='Dynamic,Tilted')
    f,p = fourier.meanPSD(figtemptilt,win=np.hanning,dx=dx,irregular=True,\
                          minpx=200)
    plt.loglog(f,p/f[0],label='Temporal,Tilted')
    f,p = fourier.meanPSD(figtempnull,win=np.hanning,dx=dx,irregular=True,\
                          minpx=200)
    plt.loglog(f,p/f[0],label='Temporal,Nulled')
    plt.legend(loc='lower left')
    plt.title('Solar B Figure - 32 Averages')
    plt.xlabel('Frequency (1/mm)')
    plt.ylabel('Power ($\mu$m$^2$ mm)')
    plt.grid()

    #Make parroting repeatability plots
    fig = plt.figure('Parroting')
    fig.add_subplot(2,2,1)
    plt.imshow(repdyntilt)
    plt.title('Dynamic Repeatability')
    plt.colorbar()
    fig.add_subplot(2,2,2)
    plt.imshow(reptemptilt)
    plt.title('Temporal Repeatability')
    plt.colorbar()
    fig.add_subplot(2,2,3)
    res = legendre2d(repdyntilt,xo=3,yo=3)[0]
    plt.imshow(repdyntilt-res)
    plt.title('Dynamic Repeatability Filtered')
    plt.colorbar()
    fig.add_subplot(2,2,4)
    res = legendre2d(reptemptilt,xo=3,yo=3)[0]
    plt.imshow(reptemptilt-res)
    plt.title('Temporal Repeatability Filtered')
    plt.colorbar()
    
def slumpedRepeatability():
    """
    Make plots demonstrating the repeatability of 4D
    metrology of slumped glass optics
    Use data taken on sample 13 on 170106
    """
    ddir = '/home/rallured/Dropbox/AXRO/Metrology/PCO1S13/170106/'
    #Plot repeatability PSDs for the 32,64,128,256 averaging
    plt.figure('Slumped Repeatability')
    plt.clf()
    for N in [32,64,128,256]:
        d1,dx = met.readCylScript(ddir+'170106_PCO1S13_Figure_DynMode_%i_0' % N)
        d2,dx = met.readCylScript(ddir+'170106_PCO1S13_Figure_DynMode_%i_1' % N)
        rep = d1-d2
        f,p = fourier.meanPSD(rep,win=np.hanning,dx=dx,irregular=True,minpx=400)
        p = p/2. #Assume uncorrelated repeatability signals
        midfreq = fourier.computeFreqBand(f,p,.1,1.,.01)*1000
        plt.loglog(f,p/f[0],label=str(N)+': %.2f nm' % midfreq)
    plt.legend(loc='lower left')
    plt.grid()
    plt.title('PCO1S13 Repeatability')
    plt.xlabel('Frequency (1/mm)')
    plt.ylabel('Power ($\mu$m$^2$ mm)')
