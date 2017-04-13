import numpy as np
import matplotlib.pyplot as plt
import axro.evaluateMirrors as eva
import astropy.io.fits as pyfits
import utilities.imaging.man as man
import legendremod as leg
import utilities.fourier as fourier
import pdb

#Load X-ray test IFs
ifs = pyfits.getdata('/home/rallured/Dropbox/AXRO/'
                     'InfluenceFunctions/220mmRoC/Cylindrical/'
                     'HIGHFIDELITY/170201_HFDFC_HighFidelity.fits')
#Get uniform stress IF
uniIF0 = pyfits.getdata('/home/rallured/Dropbox/AXRO/InfluenceFunctions/'
                       '220mmRoC/Cylindrical/Strips/170202_ConicalUniform.fits')

#Get set of IFs with uniform stress included
ifs2 = np.zeros((np.shape(ifs)[0]+1,np.shape(ifs)[1],np.shape(ifs)[2]))
ifs2[:-1] = ifs
ifs2[-1] = -uniIF0

#Get strips
strips = pyfits.getdata('/home/rallured/Dropbox/AXRO/InfluenceFunctions/220mmRoC/'
                        'Cylindrical/Strips/170202_Strips.fits')

def shademask(azwidth,axwidth):
    """
    Create a shademask over the IF shape with
    given axial and azimuthal aperture widths
    """
    shade = np.ones(np.shape(ifs[0]))
    ind = int(round(100-azwidth))
    if ind>0:
        shade[:,:ind] = 0.
        shade[:,-ind:] = 0.
    ind = int(round(100-axwidth))
    if ind>0:
        shade[:ind] = 0.
        shade[-ind:] = 0.
    return shade

def addUniform(distortion,coeff):
    """
    Rescale the uniform stress IF and add to
    distortion image
    """
    uniIF = man.newGridSize(uniIF0,np.shape(distortion),method='cubic')
    return distortion+coeff*uniIF

def correctFigurePlusSag(distortion,dx,sag=.2424,uniform=0.,\
                         azwidth=20.,axwidth=100.,smax=5.,tifs=ifs,\
                         bounds=None):
    """
    Compute the corrected figure of a mirror with uniform sag error.
    Return the corrected figure
    """
    #Create distortion map
    x,y = man.autoGrid(distortion)
    sag = -leg.singleorder(x,y,0,2)/1.5*sag #0.25 micron sag
    distortion = distortion+sag
    uniIF = man.newGridSize(uniIF0,np.shape(distortion),method='cubic')
    distortion = distortion+uniform*uniIF
    print 'Distortion ok'

    #Create shademask
    shade = shademask(azwidth,axwidth)

    #Run solver
    correction,volt = eva.correctXrayTestMirror(distortion,tifs,\
                                           shade=shade,dx=[dx],smax=smax,\
                                                bounds=bounds)
    distortion[np.isnan(correction)] = np.nan
    print 'Correction ok'

    #Run performance predictions
    pre = eva.computeMeritFunctions(distortion,[dx],renorm=True)
    post = eva.computeMeritFunctions(distortion+correction,[dx],renorm=True)
    print 'Evaluation ok'

    return distortion,correction,pre,post,volt

def uniformStressWithMask(distortion,mask,dx,uniform=0.,azwidth=70.,\
                       axwidth=70.,smax=5.):
    """
    Scan through stress values and determine corrected HPD
    vs. stress
    """
    #Add stress
    uniIF = man.newGridSize(uniIF0,np.shape(distortion),method='cubic')
    distortion2 = distortion+uniform*uniIF
    mask2 = mask+uniform*uniIF

    #Run correction
    shade = shademask(azwidth,axwidth)
    correction,volt = eva.correctXrayTestMirror(mask2,ifs,\
                                shade=shade,dx=[dx],smax=smax)
    print volt[-1]

    #Compute performance
    perf = eva.computeMeritFunctions(distortion2+correction,[dx],renorm=True)[1]

    return perf

def plotResults(res,fig=None):
    """
    Plot up results of mirror correction.
    """
    if fig is None:
        fig = plt.figure()
    fig.add_subplot(1,3,1)
    p = plt.imshow(man.remove2DLeg(man.stripnans(res[0]),xo=10))
    p.axes.set_yticks([])
    p.axes.set_xticks([])
    plt.title('PCO1S08 Figure')
    plt.colorbar()
    fig.add_subplot(1,3,2)
    p = plt.imshow(man.remove2DLeg(man.stripnans(res[1]),xo=10))
    p.axes.set_yticks([])
    p.axes.set_xticks([])
    plt.title('Flattened Piezo Correction')
    plt.colorbar()
    fig.add_subplot(1,3,3)
    p = plt.imshow(man.remove2DLeg(man.stripnans(res[0]+res[1]),xo=10))
    p.axes.set_yticks([])
    p.axes.set_xticks([])
    plt.title('Flattened Corrected Surface')
    plt.colorbar()

def fullAnalysis(img,mask,dx,azwidth=20.,axwidth=100.,res=None,filebase=None):
    """
    1) Raw correctability with 15 mm mask and 25 mm mask
    2) Masked correctability with 15 mm mask and 25 mm mask
    3) Raw and masked correctability without upper voltage bound
    4) Raw and masked aperture PSF baseline
    """
    if res is None:
        #PSF baseline
        baseline = np.copy(img)
        baseline[~np.isnan(img)] = 0.
        shade = shademask(azwidth,axwidth)
        baseline[shade==0] = np.nan
        baselineHPD = eva.computeMeritFunctions(baseline,[dx])[1]

        baseline[np.isnan(mask)] = np.nan
        maskedbaselineHPD = eva.computeMeritFunctions(baseline,[dx])[1]

        #Raw correctability
        raw = correctFigurePlusSag(img,dx,sag=0.,\
                                   azwidth=azwidth,axwidth=axwidth)
        masked = correctFigurePlusSag(mask,dx,sag=0.,\
                                      azwidth=azwidth,axwidth=axwidth)
        
        #Masked performance with dimples
        maskedwithdimple = eva.computeMeritFunctions(masked[1]+raw[0],[dx])[1]

        #No voltage bound
        raw2 = correctFigurePlusSag(img,dx,sag=0.,\
                                   azwidth=azwidth,axwidth=axwidth,smax=100.)
        masked2 = correctFigurePlusSag(mask,dx,sag=0.,\
                                      azwidth=azwidth,axwidth=axwidth,smax=100.)

        #Masked performance with dimples
        masked2withdimple = eva.computeMeritFunctions(raw2[0]+masked2[1],[dx])[1]

    else:
        [baselineHPD,maskedbaselineHPD,raw,masked,raw2,masked2,\
            maskedwithdimple,masked2withdimple] = res

    try:
        #Make plots
        fsize = (18,4)
        fig = plt.figure('raw',figsize=fsize)
        fig.add_subplot(1,3,1)
        p = plt.imshow(man.remove2DLeg(man.stripnans(raw[0]),xo=10))
        p.axes.set_yticks([])
        p.axes.set_xticks([])
        plt.title('Figure')
        plt.colorbar()
        fig.add_subplot(1,3,2)
        p = plt.imshow(man.remove2DLeg(man.stripnans(raw[1]),xo=10))
        p.axes.set_yticks([])
        p.axes.set_xticks([])
        plt.title('Flattened Piezo Correction')
        plt.colorbar()
        fig.add_subplot(1,3,3)
        p = plt.imshow(man.remove2DLeg(man.stripnans(raw[0]+raw[1]),xo=10))
        p.axes.set_yticks([])
        p.axes.set_xticks([])
        plt.title('Flattened Corrected Surface')
        plt.colorbar()
        plt.savefig(filebase+'_Raw.png')

        #Masked
        fig = plt.figure('masked',figsize=fsize)
        fig.add_subplot(1,3,1)
        p = plt.imshow(man.remove2DLeg(man.stripnans(raw[0]),xo=10))
        p.axes.set_yticks([])
        p.axes.set_xticks([])
        plt.title('Figure')
        plt.colorbar()
        fig.add_subplot(1,3,2)
        p = plt.imshow(man.remove2DLeg(man.stripnans(masked[1]),xo=10))
        p.axes.set_yticks([])
        p.axes.set_xticks([])
        plt.title('Flattened Piezo Correction')
        plt.colorbar()
        fig.add_subplot(1,3,3)
        p = plt.imshow(man.remove2DLeg(man.stripnans(raw[0]+masked[1]),xo=10))
        p.axes.set_yticks([])
        p.axes.set_xticks([])
        plt.title('Flattened Corrected Surface')
        plt.colorbar()
        plt.savefig(filebase+'_Masked.png')

        #Mid freq contribution
        img1 = man.stripnans(raw[1]+raw[0])
        img2 = man.remove2DLeg(img1,xo=10,yo=10)
        fig = plt.figure('midfreq',figsize=fsize)
        fig.add_subplot(1,3,1)
        f,p = fourier.meanPSD(img1)
        plt.loglog(f,p/f[0],label='Unfiltered')
        f,p = fourier.meanPSD(img2)
        plt.loglog(f,p/f[0],label='High Pass')
        f,p = fourier.meanPSD(img1-img2)
        plt.loglog(f,p/f[0],label='Low Pass')
        plt.grid()
        plt.xlabel('Frequency (1/mm)')
        plt.ylabel('Power ($\mu$m$^2$ mm)')
        plt.title('Mean PSD')
        plt.legend(loc='lower left')
        fig.add_subplot(1,3,2)
        p = plt.imshow(man.remove2DLeg(img1-img2,xo=10))
        p.axes.set_yticks([])
        p.axes.set_xticks([])
        plt.title('Filtered Correction')
        plt.colorbar()
        fig.add_subplot(1,3,3)
        p = plt.imshow(man.remove2DLeg(img1,xo=10))
        p.axes.set_yticks([])
        p.axes.set_xticks([])
        plt.title('Flattened Corrected Figure')
        plt.colorbar()
        plt.savefig(filebase+'_MidFreq.png')
        
    except:
        pass

    return [baselineHPD,maskedbaselineHPD,raw,masked,raw2,masked2,\
            maskedwithdimple,masked2withdimple]

def makeCorPlot(raw,masked,filebase):
    fig = plt.figure('masked',figsize=(18,4))
    fig.add_subplot(1,3,1)
    p = plt.imshow(man.remove2DLeg(man.stripnans(raw[0]),xo=10),vmin=-.3)
    p.axes.set_yticks([])
    p.axes.set_xticks([])
    plt.title('Figure')
    plt.colorbar()
    fig.add_subplot(1,3,2)
    p = plt.imshow(man.remove2DLeg(man.stripnans(masked[1]),xo=10))
    p.axes.set_yticks([])
    p.axes.set_xticks([])
    plt.title('Flattened Piezo Correction')
    plt.colorbar()
    fig.add_subplot(1,3,3)
    p = plt.imshow(man.remove2DLeg(man.stripnans(raw[0]+masked[1]),xo=10),vmin=-.3)
    p.axes.set_yticks([])
    p.axes.set_xticks([])
    plt.title('Flattened Corrected Surface')
    plt.colorbar()
    plt.savefig(filebase+'_Masked.png')
    return

def stressCompensation():
    """
    Determine whether azimuthally varying sag can produce the
    distortions measured in S21.
    Determine best amount of idealized sag to put into the mirror
    in order to get correctability?
    """
    #Get distortion data
    bounds = []
    for i in range(112):
        bounds.append([0.,5.])
    bounds.append([-100.,100.])

    perf = correctFigurePlusSag(np.zeros((200,200)),[.5],tifs=ifs2,sag=-30.,\
                                azwidth=70.,axwidth=70.,bounds=bounds)
    
    return perf
