import numpy as np
import matplotlib.pyplot as plt
import axro.solver as slv
import pyfits
import pdb
from zernikemod import stripnans
from utilities import fourier
import scipy.interpolate as interp
from astropy.modeling import models, fitting
import os
from utilities.imaging.analysis import rms

def flatSampleIF(filename,Nx,Ny,method='cubic'):
    """Read in CSV data from Vanessa and form a 2D array
    Interpolate onto grid of Nx and Ny points, where
    x is axial and y is azimuthal
    Axial is the 5mm cell direction
    Returns 2D array
    """
    #Read in data
    d = np.transpose(np.genfromtxt(filename,skip_header=1,delimiter=','))
    x = d[2]+d[5]
    y = d[3]+d[6]
    z = (d[4]+d[7])*1e6

    #Interpolate onto appropriate grid
    gx = np.linspace(y.min(),y.max(),Nx)
    gy = np.linspace(x.min(),x.max(),Ny)
    gx,gy = np.meshgrid(gx,gy)
    d = np.transpose(interp.griddata((x,y),z,(gy,gx),method=method))
    d[np.isnan(d)] = 0.
    
    return d

#Global directory variables for problem
datadir = '/home/rallured/data/solve_pzt/'

def createDist(amp,freq,phase,filename):
    """Creates a sinusoidal ripple distortion
    map for DFC2. Saves in usual ~/data directory
    """
    #Create distortion array
    y,x = np.mgrid[0:128,0:128]
    x,y = x*(100./128),y*(100./128)
    d = amp*np.sin(2*np.pi*freq*y+phase)

    #Save as fits file
    pyfits.writeto(filename,d,clobber=True)

    return d

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
    pyfits.writeto(filename,d,clobber=True)

    return d

def flatCorrection(amp,freq,phase,dx=100./128):
    """Use this function to investigate DFC MTF"""
    #Set distortion array
    d = createDist(amp,freq,phase,'dfcdist.fits')
    shade = pyfits.getdata(datadir+'shademasks/DFCmask2.fits')
    d2 = np.copy(d)
    d2[shade==0] = np.nan
    d2 = stripnans(d2)
    d2 = d2 - np.nanmean(d2)

    #Run solver to create solution
##    os.chdir(datadir+'parfiles')
##    os.system(solvedir+'solve_pzt @FlatFigure.par '
##              'premath=RoundMath.dat')
    distortionf = 'dfcdist.fits'#datadir+'distortions/dfcdist.fits'
    shadef = datadir+'shademasks/DFCmask2.fits'
##    ifuncf = '/home/rallured/Dropbox/WFS/SystemAlignment/DFC2/150624IFs/150804_RescaledIFs.fits'
    ifuncf = '/home/rallured/data/solve_pzt/ifuncs/FlatFigureMirror/150728_Resampled.fits'
    res = slv.slopeOptimizer2(ifuncf=ifuncf,distortionf=distortionf,\
                              shadef=shadef,dx=dx,smax=5.)


    #Load solution and ignore masked region
    resid = res[0] - d
    resid[shade==0] = np.nan
    resid = stripnans(resid)
    resid = resid - np.nanmean(resid)

    #Get windowed PSDs
    f,axpsdw = fourier.realPSD(resid,win=np.hanning,dx=dx)
    f = f[0] #Select only axial frequencies
    w = 2*np.pi*f/1000.
    axpsdw = axpsdw[:,0] #Select only axial frequencies

    f,origpsdw = fourier.realPSD(d2,win=np.hanning,dx=dx)
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

def toleranceEffect(shade=False):
    """Analyze effects of bonding misalignments on DFC IFs"""
    tx45 = flatSampleIF(datadir+'ifuncs/FlatFigureMirror/Tolerances/'
                      '5mmx1cm_IF_Act_45_2mmX.csv',150,150)
    ty45 = flatSampleIF(datadir+'ifuncs/FlatFigureMirror/Tolerances/'
                      '5mmx1cm_IF_Act_45_2mmY.csv',150,150)
    tr45 = flatSampleIF(datadir+'ifuncs/FlatFigureMirror/Tolerances/'
                      '5mmx1cm_IF_Act_45_2degCW.csv',150,150)
    tx32 = flatSampleIF(datadir+'ifuncs/FlatFigureMirror/Tolerances/'
                      '5mmx1cm_IF_Act_32_2mmX.csv',150,150)
    ty32 = flatSampleIF(datadir+'ifuncs/FlatFigureMirror/Tolerances/'
                      '5mmx1cm_IF_Act_32_2mmY.csv',150,150)
    tr32 = flatSampleIF(datadir+'ifuncs/FlatFigureMirror/Tolerances/'
                      '5mmx1cm_IF_Act_32_2degCW.csv',150,150)

    if shade==True:
        shade = pyfits.getdata(datadir+'shademasks/roundmask3.fits')
        tx45[shade==0] = np.nan
        ty45[shade==0] = np.nan
        tr45[shade==0] = np.nan
        tx32[shade==0] = np.nan
        ty32[shade==0] = np.nan
        tr32[shade==0] = np.nan
        

    ifunc = pyfits.getdata(datadir+'ifuncs/FlatFigureMirror/'
                           '150319FlatIFs.fits')*1e3

    #Make plots
    fig = plt.figure()
    fig.add_subplot(231)
    plt.imshow(tx45-ifunc[44,0])
    plt.colorbar()
    plt.title('2mm X Trans - Cell 45')
    fig.add_subplot(232)
    plt.imshow(ty45-ifunc[44,0])
    plt.colorbar()
    plt.title('2mm Y Trans - Cell 45')
    fig.add_subplot(233)
    plt.imshow(tr45-ifunc[44,0])
    plt.colorbar()
    plt.title('2 deg Rot - Cell 45')
    fig.add_subplot(234)
    plt.imshow(tx32-ifunc[31,0])
    plt.colorbar()
    plt.title('2mm X Trans - Cell 32')
    fig.add_subplot(235)
    plt.imshow(ty32-ifunc[31,0])
    plt.colorbar()
    plt.title('2mm Y Trans - Cell 32')
    fig.add_subplot(236)
    plt.imshow(tr32-ifunc[31,0])
    plt.colorbar()
    plt.title('2 deg Rot - Cell 32')

def interpolatedFilter(model,dmodel,measurement,dmeas):
    """Interpolate log of FEA FFT component magnitudes onto measured
    IF frequency grid. Take the ratio as a 2D frequency filter. Fit this
    to a low order 2D polynomial and use this model as a filter."""
    #Normalize model to measurement
    model = model * (measurement.max()-measurement.min())/\
            (model.max()-model.min())
    #Generate Fourier magnitudes of model and measurement
    modelFFT = np.abs(fourier.continuousComponents(model,dmodel))
    measFFT = np.abs(fourier.continuousComponents(measurement,dmeas))
    #Construct frequency grids
    modelfx,modelfy = fourier.freqgrid(model)
    measfx,measfy = fourier.freqgrid(measurement)
    #Use griddata to interpolate the log of the model onto the
    #frequency grid of the measurement
    #Linear method is good enough given that we will
    #fit the filter with a low order polynomial
    newmodel = interp.griddata((modelfx.flatten(),modelfy.flatten()),\
                        np.log10(modelFFT.flatten()),\
                        (measfx,measfy),method='linear')

    #Select positive frequency region of FFT
    #Set constant term to 0., and normalize to a maximum
    #log ratio of 0.
    lograt = newmodel - np.log10(measFFT)
    lograt[0,0] = 0.
    lograt[lograt>0.] = 0.
    sh = np.shape(lograt)
    lograt = lograt[:sh[0]/2,:sh[1]/2]
    fxsub = measfx[:sh[0]/2,:sh[1]/2]
    fysub = measfy[:sh[0]/2,:sh[1]/2]
    #Perform quadratic fit to the log ratio
    p_init = models.Polynomial2D(degree=2)
    fit_p = fitting.LevMarLSQFitter()
    p = fit_p(p_init,fxsub,fysub,lograt)
    #Form a function covering full frequency range
    lograt = p(np.abs(measfx),np.abs(measfy))
    ratio = 10**lograt
    
    return measFFT,ratio

def spiePlot():
    """Make up ripple induction plot for SPIE 2015 slides"""
    #Load data
    os.chdir('/Users/ryanallured/Dropbox/WFS/'
             'SystemAlignment/DFC2/Iteration_0_3/ForAlexey')
    dist = pyfits.getdata('dfcdist_0_3.fits')
    shade = pyfits.getdata('Shademask.fits')
    pred = np.genfromtxt('Predict.txt')
    p2 = np.genfromtxt('Phase02.txt')

    #Apply shademask
    dist[shade==0] = np.nan
    pred[shade==0] = np.nan
    p2[shade==0] = np.nan

    #Create residual arrays
    preresid = dist-pred
    achresid = dist-p2
    indresid = pred-p2

    #Convert to axial slope
    dist = np.diff(dist,axis=0)/1000./(100./124)*180/np.pi*60**2
    preresid = np.diff(preresid,axis=0)/1000./(100./124)*180/np.pi*60**2
    achresid = np.diff(achresid,axis=0)/1000./(100./124)*180/np.pi*60**2
    indresid = np.diff(indresid,axis=0)/1000./(100./124)*180/np.pi*60**2

    #Remove nans
    dist = stripnans(dist)
    preresid = stripnans(preresid)
    achresid = stripnans(achresid)
    indresid = stripnans(indresid)

    #Make plots
    plt.figure()
    plt.subplot(141)
    plt.imshow(dist)
    plt.title('Desired')
    plt.subplot(142)
    plt.imshow(preresid,vmax=dist.max(),vmin=dist.min())
    plt.title('Predicted')
    plt.subplot(143)
    plt.imshow(achresid,vmax=dist.max(),vmin=dist.min())
    plt.title('Achieved')
    plt.subplot(144)
    plt.imshow(indresid,vmax=dist.max(),vmin=dist.min())
    plt.title('Predicted-Achieved')

    #Print out RMS axial slopes
    print rms(dist)
    print rms(preresid)
    print rms(achresid)
    print rms(indresid)

    #Look at PSD of residual
    f,p = fourier.realPSD(indresid,dx=100./124)
    f = f[0]
    plt.figure()
    plt.semilogy(f,p[:,0])
    return indresid
