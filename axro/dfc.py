import numpy as np
import matplotlib.pyplot as plt
import axro.solver as slv
import pyfits
import pdb
from zernikemod import stripnans
from utilities import fourier

#Global directory variables for problem
datadir = '/home/rallured/data/solve_pzt/'

def createDist(amp,freq,phase,filename):
    """Creates a sinusoidal ripple distortion
    map for DFC2. Saves in usual ~/data directory
    """
    #Create distortion array
    y,x = np.mgrid[0:128,0:128]
    x,y = x*(100./124),y*(100./124)
    d = amp*np.sin(2*np.pi*freq*y+phase)

    #Save as fits file
    pyfits.writeto(datadir+'distortions/'+filename,d,clobber=True)

def runSolver(ifuncf,distortionf,shadef,vfile,pfile,iterate=None):
    """Runs the Python solver and formulates the
    proper voltage array for the IRIS controller.
    This array is saved in Dropbox for application
    to the mirror"""
    #If this is an iteration, need new bounds
    if iterate != None:
        iterate = iterate/5.
        bounds = []
        for i in range(np.size(iterate)):
            if i != 2 and i != 68:
                bounds.append((-iterate[i],1.-iterate[i]))
        pdb.set_trace()
        res = slv.slopeOptimizer2(ifuncf=ifuncf,shadef=shadef,\
                            distortionf=distortionf,dx=100./120,bounds=bounds)
    else:
        #Runs the solver and gets the voltages
        res = slv.slopeOptimizer2(ifuncf=ifuncf,shadef=shadef,\
                                  distortionf=distortionf,dx=100./120,smax=1.)
    voltages = res[1] * 5.
    voltages[voltages>5.] = 5.
    
    #Insert missing voltages for cell 3 and 69
    voltages = np.insert(voltages,2,0.)
    voltages = np.insert(voltages,68,0.)
    #Writes the voltage array to the Dropbox
    np.savetxt('/home/rallured/Dropbox/WFS/SystemAlignment/DFC2/'+vfile,\
                   voltages)
    #Writes the predicted array to the Dropbox
    np.savetxt('/home/rallured/Dropbox/WFS/SystemAlignment/DFC2/'+pfile,\
                   res[0])

    return res

def createDistSim(amp,freq,phase,filename):
    """Creates a sinusoidal ripple distortion
    map for DFC2. Saves in usual ~/data directory
    """
    #Create distortion array
    y,x = np.mgrid[0:150,0:150]
    x,y = x*(100./150),y*(100./150)
    d = amp*np.sin(2*np.pi*freq*y+phase)

    #Save as fits file
    pyfits.writeto(datadir+'distortions/'+filename,d,clobber=True)

    return d

def flatCorrection(amp,freq,phase):
    """Use this function to investigate DFC MTF"""
    #Set distortion array
    d = createDistSim(amp,freq,phase,'dfcdist.fits')
    shade = pyfits.getdata(datadir+'shademasks/roundmask3.fits')
    d2 = np.copy(d)
    d2[shade==0] = np.nan
    d2 = stripnans(d2)
    d2 = d2 - np.nanmean(d2)

    #Run solver to create solution
##    os.chdir(datadir+'parfiles')
##    os.system(solvedir+'solve_pzt @FlatFigure.par '
##              'premath=RoundMath.dat')
    distortionf = datadir+'distortions/dfcdist.fits'
    shadef = datadir+'shademasks/roundmask3.fits'
    ifuncf = datadir+'ifuncs/FlatFigureMirror/150715_TwoShorted.fits'
    res = slv.slopeOptimizer2(ifuncf=ifuncf,distortionf=distortionf,\
                              shadef=shadef,dx=100./150)


    #Load solution and ignore masked region
    resid = res[0] - d
    resid[shade==0] = np.nan
    resid = stripnans(resid)
    resid = resid - np.nanmean(resid)

    #Get windowed PSDs
    f,axpsdw = fourier.realPSD(resid,win=np.hanning,dx=100./150)
    f = f[0] #Select only axial frequencies
    w = 2*np.pi*f/1000.
    axpsdw = axpsdw[:,0] #Select only axial frequencies

    f,origpsdw = fourier.realPSD(d2,win=np.hanning,dx=100./150)
    f = f[0]
    origpsdw = origpsdw[:,0]
    
##    ripple = sum((w**2*axpsdw)[f>.15])/sum(w**2*origpsdw)

##    #Get unwindowed PSDs for input reduction
##    f,axpsd = fourier.realPSD(resid,dx=.5)
##    f = f[0] #Select only axial frequencies
##    axpsd = axpsd[:,0] #Select only axial frequencies
##
##    f,origpsd = fourier.realPSD(d,dx=.5)
##    f = f[0]
##    origpsd = origpsd[:,0]
    
##    correction = sum((w**2*axpsdw)[f<.15])/sum(w**2*origpsdw)
    correction = sum(w**2*axpsdw)/sum(w**2*origpsdw)
    
    return correction
