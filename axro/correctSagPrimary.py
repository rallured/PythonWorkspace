import numpy as np
import matplotlib.pyplot as plt
import traces.axro.slf as slf
import axro.evaluateMirrors as eva
import pyfits
import utilities.imaging.man as man
import utilities.imaging.fitting as fit
import legendremod as leg
import pdb
import scipy.signal as sig

#Load IFs
ifs = pyfits.getdata('/home/rallured/Dropbox/AXRO/'
                     'XrayTest/IFs/161019_5x5mm_IFs.fits')

def correctSag(azwidth=20.,axwidth=100.):
    """
    Compute the corrected figure of a mirror with uniform sag error.
    Convert the corrected figure into 2D Legendre coefficients.
    These coefficiens can be used as input to the SLF raytrace.
    """
    #Create distortion map
    x,y = man.autoGrid(ifs[0])
    sag = -leg.singleorder(x,y,0,2)/1.5*.2424 #0.25 micron sag

    #Create shademask
    shade = np.ones(np.shape(ifs[0]))
    ind = int(round(100-azwidth))
    if ind>0:
        shade[:,:ind] = 0.
        shade[:,-ind:] = 0.
    ind = int(round(100-axwidth))
    if ind>0:
        shade[:ind] = 0.
        shade[-ind:] = 0.

    #Run solver
    correction = eva.correctXrayTestMirror(sag,ifs,shade=shade,dx=[.5,.5])[0]

    #Fit Legendres to result
    res = man.stripnans(sag+correction)
    coeff,ax,az = fit.fitLegendreDistortions(res,xo=20,yo=20)

    #coefficients in mm, axial order, azimuthal order
    #Get rid of anything 5 nm or less in coefficient
    ind = np.abs(coeff)>.005
    coeff = coeff[ind]
    ax = ax[ind]
    az = az[ind]

    return [coeff/1000.,ax,az],res

def determinePerformance(plist=[[0],[0],[0]],azwidth=20.,axwidth=100.):
    """
    Use the coefficients from correctSag as input to the
    SLF raytrace. Determine optimal performance using
    pitch compensation.
    """
    #Run traces for various pitch compensation
    pitch = np.linspace(0.,10*.3e-3,200)
    res = np.transpose(np.array([slf.singleOptic2(1000,misalign=[0,0,0,0,t,0],\
                            az=azwidth,ax=axwidth,\
                            plist=plist) for t in pitch]))
    #Fit optimal performance
    per = sig.savgol_filter(res[0],11,3)
    ind = np.argmin(per)
    return [res[0][ind],pitch[ind],res[1][ind]]
