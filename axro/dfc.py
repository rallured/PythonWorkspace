import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import axro.solver as slv
import pyfits
import pdb
from zernikemod import stripnans
from utilities import fourier
import scipy.interpolate as interp
from astropy.modeling import models, fitting
import os
from utilities.imaging.analysis import rms
import utilities.imaging.man as man
import utilities.imaging.fitting as ufit
from reconstruct import reconstruct
import utilities.imaging.analysis as anal
import utilities.imaging.fitting as fitting
from scipy.interpolate import griddata

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
    #d[np.isnan(d)] = 0.
    
    return d

def simWFSSmooth(ifs):
    """Simulate the resolution effect of the WFS. Take the
    derivatives of the IFs, rebin to the WFS lenslet size,
    integrate the slopes into a reconstructed wavefront with
    the Southwell algorithm, then run these new IFs through
    the CTF code
    """
    #Compute slopes, rebin to 128 x 128 pixels, integrate
    #to get new IF
    N = np.shape(ifs)[0]
    newifs = np.zeros((N,128,128))
    for i in range(N):
        #Get slopes
        xang = np.diff(ifs[i],axis=0)/(100./513*1000)
        xang = xang[:,:512]
        yang = np.diff(ifs[i],axis=1)/(100./513*1000)
        yang = yang[:512,:]
        #Logical or NaN indices
        ind = np.logical_or(np.isnan(xang),np.isnan(yang))
        xang[ind] = np.nan
        yang[ind] = np.nan
        #Rebin appropriately
        xang = man.rebin(xang,(128,128))
        yang = man.rebin(yang,(128,128))
        #Pad with NaNs
        xang = man.padRect(xang)
        yang = man.padRect(yang)
        #Create phase array
        ph = np.zeros(np.shape(xang))
        ph[np.isnan(xang)] = 100.
        xang[np.isnan(xang)] = 100.
        yang[np.isnan(yang)] = 100.
        #Format for Fortran
        ph = np.array(ph,order='F')
        xang = np.array(xang,order='F')
        yang = np.array(yang,order='F')
        #Reconstruct into phase
        ph2 = reconstruct(xang,yang,1e-12,100./128*1000,ph)
        ph2[ph2==100] = np.nan
        newifs[i] = stripnans(ph2)

    return newifs

    

#Global directory variables for problem
datadir = '/home/rallured/data/solve_pzt/'

def createDist(amp,freq,phase,filename):
    """Creates a sinusoidal ripple distortion
    map for DFC2. Saves in usual ~/data directory
    """
    #Create distortion array
    y,x = np.mgrid[0:128,0:128]
    x,y = x*(100./125),y*(100./125)
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
    ifuncf = '/home/rallured/Dropbox/WFS/SystemAlignment/DFC2/150730IFs/150917_Gauss5.fits'#150917_Gauss2.fits'
##    ifuncf = '/home/rallured/data/solve_pzt/ifuncs/FlatFigureMirror/150728_Resampled.fits'
    res = slv.slopeOptimizer2(ifuncf=ifuncf,distortionf=distortionf,\
                              shadef=shadef,dx=dx,smax=100.)


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

    #Get unwindowed PSDs for input reduction
    f,axpsd = fourier.realPSD(resid,dx=dx)
    f = f[0] #Select only axial frequencies
    axpsd = axpsd[:,0] #Select only axial frequencies

    f,origpsd = fourier.realPSD(d2,dx=dx)
    f = f[0]
    origpsd = origpsd[:,0]
    
##    correction = sum((w**2*axpsdw)[f<.15])/sum(w**2*origpsdw)
##    correction = sum(w**2*axpsdw)/sum(w**2*origpsdw)
    correction = (anal.rms(np.diff(resid,axis=0))/anal.rms(np.diff(d2,axis=0)))**2
    
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

#Examine peak locations
def peakLocations(fileIF):
    """This examines the IFs measured on DFC2 to determine
    if there is an offset in the peak locations. This could
    potentially explain the low frequency distortions we are
    seeing in the corrections."""
    #Load in influence functions
    d = pyfits.getdata(fileIF)

    #Loop through and get indices of max pixel
    N = np.shape(d)[0]
    cx = np.zeros(N)
    cy = np.zeros(N)
    for i in range(N):
        cy[i],cx[i] = np.where(d[i]==np.nanmax(d[i]))

    #Plot up the peak locations
    fig = plt.figure()
    plt.plot(cx,cy,'*')

#Examine locality of measured vs. FEA IFs
def locality():
    """Examine the locality of the IFs for DFC2. IFs should
    be in (N,128,128) shape, where first index matches cell
    to cell for each IF.
    Power spectra are computed for each cell, with nearest
    neighbor interpolation for the measured IFs.
    The PSDs are then normalized to sum to unity and then
    the expected value for frequency taking the PSD to be
    a distribution is taken as the locality figure of merit.
    A large figure of merit indicates a well-localized IF.
    Return FoM ratios for each cell.
    """
    #Load in the IFs
    fea = pyfits.getdata('/home/rallured/data/solve_pzt/ifuncs/'
                         'FlatFigureMirror/150916_128binswithNaNs.fits')
    meas = pyfits.getdata('/home/rallured/Dropbox/WFS/SystemAlignment/'
                          'DFC2/150730IFs/150914_RescaledIFs.fits')
    N = np.shape(meas)[0]

    #Fill in NaNs
    #Loop through, compute PSDs, compute widths
    measwidth = []
    feawidth = []
    for i in range(N):
##        meas[i] = man.nearestNaN(meas[i])
##        measf,measp = fourier.realPSD(meas[i],dx=100./125)
##        measp = np.sqrt(measp)
##        measp = measp/np.sum(measp)
##        fx,fy = np.meshgrid(measf[1],measf[0])
##        fr = np.sqrt(fx**2+fy**2)
##        fr = fr.flatten()
##        measp = measp.flatten()
##        measwidth.append(np.average(fr,weights=measp))
##        feaf,feap = fourier.realPSD(fea[i],dx=100./128)
##        feap = np.sqrt(feap)
##        feap = feap/np.sum(feap)
##        fx,fy = np.meshgrid(feaf[1],feaf[0])
##        fr = np.sqrt(fx**2+fy**2)
##        fr = fr.flatten()
##        feap = feap.flatten()
##        feawidth.append(np.average(fr,weights=feap))
        #Compute using image moments
        fea[i] = fea[i] - np.median(fea[i][np.invert(np.isnan(fea[i]))])
        cx,cy,stdx,stdy = anal.findMoments(fea[i]**2)
        feawidth.append(np.sqrt(stdx**2+stdy**2))
        meas[i] = meas[i] - np.median(meas[i][np.invert(np.isnan(meas[i]))])
        cx,cy,stdx,stdy = anal.findMoments(meas[i]**2)
        measwidth.append(np.sqrt(stdx**2+stdy**2))
        pdb.set_trace()

    return np.array(measwidth),np.array(feawidth)
    
def localizationEx():
    #Load in the IFs
    fea = pyfits.getdata('/home/rallured/data/solve_pzt/ifuncs/'
                         'FlatFigureMirror/150728_Resampled.fits')
    meas = pyfits.getdata('/home/rallured/Dropbox/WFS/SystemAlignment/'
                          'DFC2/150730IFs/150914_RescaledIFs.fits')
    meas[43] = man.nearestNaN(meas[43])

    #PSDs
    fm,pm = fourier.realPSD(meas[43],dx=100./125)
    ff,pf = fourier.realPSD(fea[43],dx=100./128)

    #Make plot
    fig = plt.figure()
    fig.add_subplot(121)
    plt.imshow(pf,interpolation='none',norm=LogNorm(),\
               extent=[ff[0][0],ff[0][-1],ff[0][-1],ff[0][0]])
    plt.title('Modeled IF - Cell 45')
    plt.xlabel('Azimuthal Frequency (1/mm)')
    plt.ylabel('Axial Frequency (1/mm)')
    plt.colorbar()
    fig.add_subplot(122)
    plt.imshow(pm,interpolation='none',norm=LogNorm(),\
               extent=[fm[0][0],fm[0][-1],fm[0][-1],fm[0][0]])
    plt.title('Measured IF - Cell 45')
    plt.xlabel('Azimuthal Frequency (1/mm)')
    plt.ylabel('Axial Frequency (1/mm)')
    plt.colorbar()

    fig = plt.figure()
    sl = meas[43][:,70]
    fsl = fea[43][:,70]
    fsl = fsl*(np.nanmax(sl)-np.nanmin(sl))/(np.nanmax(fsl)-np.nanmin(fsl))
    fsl = fsl - np.nanmax(fsl) + np.nanmax(sl)
    plt.plot(sl,label='Measured')
    plt.plot(np.arange(128)+3,fsl,label='Modeled')
    plt.legend(loc='upper right')
    plt.title('Central Slice of Cell 45')

    fig = plt.figure()
    plt.plot(np.diff(sl))
    plt.plot(np.arange(127)+3,np.diff(fsl))

#Need to compare repeatability and ultimate correction PSDs
#Are we matching low frequency repeatability limit? Or is there
#residual low frequency distortions that we are unable to correct?
#If there are, then can these be explained with the CTF?
def compareResidual(resid,rep,shade,dx=100./125.*1000.,win=1):
    """Create averaged slope PSDs of both the residual
    and the repeatability data. Do this both ways: using diff
    and using w**2 weighting"""
    resid[shade==0] = np.nan
    resid = man.stripnans(resid)
    resd = np.diff(resid,axis=0)
    resd = resd-np.nanmean(resd)
    
    rep[shade==0] = np.nan
    rep = man.stripnans(rep)
    repd = np.diff(rep,axis=0)
    repd = repd-np.nanmean(repd)

    #Create PSDs
    f,p = fourier.meanPSD(resid,dx=dx,win=win)
    fsl,psl = fourier.meanPSD(resd/dx,\
                              dx=dx,win=win)
    frep,prep = fourier.meanPSD(rep,dx=dx,win=win)
    fslrep,pslrep = fourier.meanPSD(repd/dx,\
                              dx=dx,win=win)
    p = p * (2*np.pi*f)**2
    prep = prep * (2*np.pi*frep)**2

    pdb.set_trace()

    #Make plots
    plt.figure()
    plt.loglog(fsl,psl/fsl[1],label='Resid, BF')
    print np.sqrt(np.sum(psl))
    print np.sqrt(np.sum(p))
    print np.sqrt(np.sum(pslrep))
    print np.sqrt(np.sum(prep))
    plt.loglog(f,p/f[1],label='Resid, W')
    plt.loglog(fslrep,pslrep/fslrep[1],label='Rep, BF')
    plt.loglog(frep,prep/frep[1],label='Rep, W')

    return [f,p],[fsl,psl],[frep,prep],[fslrep,pslrep]

def applySG(ifs,shade2,n,m):
    """Apply a Savitzky-Golay filter to a set of IFs
    Shade2 should be a shade mask slightly bigger than
    the intended shade mask in order to contain all
    needed data rows when taking the slopes
    n is the window size (odd)
    and m is the filter order
    """
    #Apply shade mask and SG filter
    sh = np.shape(ifs)
    for i in range(1,sh[0]):
        ifs[i][shade2==0] = np.nan
        ifs[i][shade2==1] = fitting.sgolay2d(\
            man.stripnans(ifs[i]),n,m)[0].flatten()

    return ifs

def SGexperiment(d,dx,shade2,n,m,noise,ndx):
    """Conduct an experiment to determine optimal SG filter.
    Need a theoretical (no noise) IF and a noise model (metrology data).
    SG filter cannot handle NaNs, so shade mask is required.
    Add noise model to data, Apply SG filter, compare result to original
    IF. Continue applying SG filter and determine how many repetitions
    produce the best result.
    Just take one instance of repeatability data as noise
    Metric should be RMS deviation from result and original.
    This can then be repeated for varying n,m to determine the optimal
    SG filter parameters."""
    #Apply shademask to input data
    d = np.copy(d)
    noise = np.copy(noise)
    d[shade2==0] = np.nan
    d = man.stripnans(d)
    noise[shade2==0] = np.nan
    noise = man.stripnans(noise)

    #Load in noise model
    nfreqx,nfreqy = fourier.freqgrid(noise,dx=ndx)
    noisef = fourier.continuousComponents(noise,ndx)

    #Interpolate noise model to data frequency grid
    newx,newy = fourier.freqgrid(d,dx=dx)
    noisef2 = griddata((nfreqx.flatten(),nfreqy.flatten()),noisef.flatten(),\
                     (newx,newy),\
                     method='linear',fill_value=0.)
    #Add noise to input
    noisef2 = noisef2/dx
    newnoise = np.real(np.fft.ifftn(noisef2))
    d2 = d + newnoise

    #Apply SG Filter
    df = fitting.sgolay2d(d2,n,m)[0]

    pdb.set_trace()

    #Compare with original
    print anal.rms(df-d)
    
    pdb.set_trace()
    
    return None
