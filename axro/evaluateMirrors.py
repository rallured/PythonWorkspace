import numpy as np
import matplotlib.pyplot as plt
import scattering as scat
import utilities.imaging.man as man
import utilities.imaging.fitting as fit
import utilities.imaging.analysis as anal
import scipy.ndimage as nd
import pdb
import axro.solver as slv
from utilities.plotting import myhist
import traces.conicsolve as conic
import legendremod as leg
from utilities.imaging.analysis import readCyl4D

def correctXrayTestMirror(d,ifs,shade=None,dx=None,azweight=.015,smax=5.):
    """
    Get distortion on same grid as IFs and run correction.
    Rebin result onto original distortion grid and apply.
    dx should be on IF grid size
    """
    #Rebin to IF grid
    d2 = man.newGridSize(d,np.shape(ifs[0]))

    #Handle shademask
    if shade is None:
        shade = np.ones(np.shape(d2))

    #Run correction
    orig, cor, volt = slv.correctDistortion(d2,ifs,shade,\
                                            dx=dx,azweight=azweight,\
                                            smax=smax)
    
    #Add correction to original data
    ifs2 = ifs.transpose(1,2,0)
    cor2 = np.dot(ifs2,volt)
    cor3 = man.newGridSize(cor2,np.shape(d),method='cubic')
    #Handle shademask
    cor2[shade==0] = np.nan
    cornan = man.newGridSize(cor2,np.shape(d),method='linear')
    cor3[np.isnan(cornan)] = np.nan

    return cor3,volt

def computeMeritFunctions(d,dx,x0=np.linspace(-5.,5.,1000),\
                          graze=conic.woltparam(220.,8400.)[0]):
    """
    RMS axial slope
    Axial sag
    
    """
    #Remove NaNs
    d = man.stripnans(d)
    
    #Compute axial sag as a function of azimuth
##    sag = [anal.fitSag(t)*np.shape(d)[0]**2 for t in np.transpose(d)]

##    #Compute RMS axial slope
##    grad = np.gradient(d,*dx)
##    rmsy = anal.rms(grad[0])/1e3*180/np.pi*60**2
##    rmsx = anal.rms(grad[1])/1e3*180/np.pi*60**2
    
    #Compute PSF
    primfoc = conic.primfocus(220.,8400.)
    dx2 = x0[1]-x0[0]
##    z = np.linspace(8500.,8400.,np.shape(d)[0]) #Proper orientation for 4D
    z = np.arange(0,-np.shape(d)[0],-1)*dx[0]*np.cos(graze)+8400.
    res = np.array([scat.primaryPSF(z=z,perturb=-di/1e3,x0=x0) \
                    for di in np.transpose(d)])
    resa = np.mean(res,axis=0)
    if np.sum(resa)*dx2 < .95:
        print 'Possible sampling problem'
        print str(np.sum(resa)*dx2)

    #Compute PSF merit functions
    rmsPSF = np.sqrt(np.sum(resa*x0**2)*dx2-(np.sum(resa*x0)*dx2)**2)
##    h = np.array([np.sum(resa[np.abs(x0)<m])*dx2 for m in np.abs(x0)])
##    hpdPSF = np.abs(x0[np.argmin(np.abs(h-.5))])*2
    cdf = np.cumsum(resa)*dx2
    hpdPSF = x0[np.argmin(np.abs(cdf-.75))]-\
             x0[np.argmin(np.abs(cdf-.25))]
    

    #Make plots
##    plt.figure()
##    plt.clf()
##    slps = man.nanflatten(grad[0]/1e3*180/np.pi*60**2*2)
##    h = np.array(myhist(slps,bins=100))
##    plt.plot(h[1],h[0]/np.max(h[0]),label='Slopes')
##    plt.plot(x0/primfoc*180/np.pi*60**2,resa/np.max(resa),label='Kirchhoff')

##    slps = np.abs(slps)
##    slps.sort()
##    hpdsl = slps[len(slps)/2]*2

##    return slps,resa,\
##           sag,rmsy*2,rmsx*2,\
##           rmsPSF/primfoc*180/np.pi*60**2*2,hpdPSF/primfoc*180/np.pi*60**2,\
##           rmsy*2*2,hpdsl
    return rmsPSF/primfoc*180/np.pi*60**2*2,hpdPSF/primfoc*180/np.pi*60**2
    
def evaluate2DLeg(xo,yo,shade=np.ones((200,200)),ifs=None,coeff=1.,smax=5.):
    """
    Compute a 2D Legendre's and compute before and after
    sensitivity.
    """
    #Load IFs
    if ifs is None:
        ifs = pyfits.getdata('/home/rallured/Dropbox/AXRO/'
                             'XrayTest/IFs/161019_5x5mm_IFs.fits')
    
    #Create distortion
    x,y = man.autoGrid(ifs[0])
    dist = leg.singleorder(x,y,xo,yo)*coeff
    
    #Create correction
    cor,volt = correctXrayTestMirror(dist,ifs,dx=[100./201],shade=shade,\
                                     smax=smax)

    #Evalute performances
    perf0 = computeMeritFunctions(dist,[.5],x0=np.linspace(-5.,5.,500))[-1]
    perf1 = computeMeritFunctions(dist+cor,[.5],x0=np.linspace(-5.,5,500))[-1]

    return perf0,perf1,volt
