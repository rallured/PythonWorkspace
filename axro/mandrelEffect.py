import numpy as np
import matplotlib.pyplot as plt
import axro.evaluateMirrors as eval
import astropy.io.fits as pyfits
import utilities.imaging.man as man
import legendremod as leg

def singleMirrorTest(shade=np.ones((200,200))):
    """
    1) Load S08 distortion data
    2) Add sag due to mandrel error
    3) Apply appropriate shademask
    4) Compute Correction
    5) Evaluate performance pre and post correction
    """
    d = pyfits.getdata('/home/rallured/Dropbox/AXRO/Metrology/PCO1S08/'
                       '161108_PCO1S08_CleanedDistortionData.fits')
    dx = 0.151423438691
    ifs = pyfits.getdata('/home/rallured/Dropbox/AXRO/InfluenceFunctions/'
                         '100x100_5x5_220/161019_5x5mm_IFs.fits')
    x,y = man.autoGrid(d)
    sag = leg.singleorder(x,y,0,2)/1.5*-0.25 # 0.25 micron axial sag
    d = d+sag
    cor = eval.correctXrayTestMirror(d,ifs,shade=shade,dx=[dx])
    d[np.isnan(cor[0])] = np.nan
    dcor = d+cor[0]

    #Evaluate performance
    perf0 = eval.computeMeritFunctions(d,[dx])
    perf1 = eval.computeMeritFunctions(d+cor[0],[dx])

    return d,dcor,perf0,perf1

def plotRes(d,dcor):

    #Make plots
    fig = plt.figure('SingleMir')
    fig.clf()
    p = fig.add_subplot(1,3,1)
    plt.imshow(man.stripnans(d))
    p.axes.set_xticks([])
    p.axes.set_yticks([])
    plt.colorbar()
    plt.title('Initial Distortion\n')
    p = fig.add_subplot(1,3,2)
    plt.imshow(man.stripnans(dcor))
    p.axes.set_xticks([])
    p.axes.set_yticks([])
    plt.colorbar()
    plt.title('Corrected Figure\n')
    dflat = man.removeLegSlice(dcor,order=1)
    p = fig.add_subplot(1,3,3)
    plt.imshow(man.stripnans(dflat),vmax=.08,vmin=-.08)
    p.axes.set_xticks([])
    p.axes.set_yticks([])
    plt.colorbar()
    plt.title('Flattened and Rescaled\n')

    return None
