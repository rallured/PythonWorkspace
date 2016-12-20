import numpy as np
import matplotlib.pyplot as plt
import utilities.metrology as met
import utilities.fourier as fourier
import pdb

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
