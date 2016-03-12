import numpy as np
from scipy.ndimage import rotate
import zygo
import pdb,glob
import matplotlib.pyplot as plt
import utilities.imaging.analysis as anal
import utilities.imaging.man as man
from scipy.ndimage.filters import gaussian_filter1d
import utilities.imaging.fitting as fit

#Need to analyze Sydor metrology of glass wafers
#Load in metrology data, and compute slope distributions
#over 95 x 75 mm rectangle about the center of the
#image. Scan through a rotation angle. NaNs filled in
#with mean of image. Should be small effect on slope
#statistics. Image shape will be changed by rotation, so
#choose center of image independently of shape.

def slopeStatistics(filename):
    #Load image
    i,phase,lat = zygo.readzygo(filename)
    lat = lat*1000.
    #Fill missing values with mean value of image
    phase[np.isnan(phase)] = np.nanmean(phase)

    #Loop through rotation angles and compute
    #RMS slope over region of interest
    fom = []
    for ang in np.linspace(0.,180.,180):
        #Rotate image by appropriate angle
        d = rotate(phase,ang)
        #Transform to slope in dispersion direction
        sl = np.diff(d*1000./lat*180/np.pi*60**2,axis=0)
        #Select region of interest
        sh = np.shape(sl)
        sl = sl[round(sh[0]/2-75./2/lat):round(sh[0]/2+75./2/lat),\
                round(sh[1]/2-95./2/lat):round(sh[1]/2+95./2/lat)]
        #Compute rms of slope
        sl = np.abs(sl - np.mean(sl)) #Subtract average tilt
        sl = sl.flatten()
        sl = np.sort(sl)
        #Histogram slopes
##        y,x = np.histogram(sl,bins=np.linspace(-200.,200.,1000))
##        pdb.set_trace()
        #fom.append(anal.fwhm(x[1:],y))
        fom.append(sl[round(.875*np.size(sl))]-sl[round(.125*np.size(sl))])
    return np.array(fom)

#Loop through Sydor ascii files and plot FoM
#for each one
def sydorPlot(f):
    #Get filenames
##    f = glob.glob('/home/rallured/Dropbox/Arcus/Sydor/2015_07_22_150mm_flat_wafers/*.txt')
    #Loop through and plot FoM vs. rotation
    ang = np.linspace(0.,180.,180)
    plt.ion()
    plt.figure()
    largest = np.zeros(np.size(f))
    for i in range(np.size(f)):
        fom = slopeStatistics(f[i])
        plt.plot(ang,fom,label=f[i].split('.')[0])
        #largest[i] = np.max(fom[fom.argmin()-5:fom.argmin()+6])
    #Plot requirement
    plt.plot([0.,180],[17.7,17.7],'k--')
    plt.title('Gaussian FWHM Equivalent Widths')
    plt.xlabel('Wafer Rotation Angle')
    plt.ylabel('Width (arcsec)')
    return None#largest

def estimatePerf(f):
    """Estimate performance of active area of first set of
    nanoimprint gratings on Sydor glass
    """
    #Load image
    i,phase,lat = zygo.readzygo(f)
    lat = lat*1000.
    phase = phase*1000.

    pdb.set_trace()

    #Choose region
    subimg = anal.getSubApp(phase)

    #Fill in NaN values using spline interpolation
    subimg = man.nearestNaN(subimg,method='cubic')

    #Determine gradient
    gradx,grady = np.gradient(subimg,lat)
    gy = grady.flatten()*np.sin(1.5*np.pi/180)

    #Bin up and convolve with Gaussian
    n,b = np.histogram(gy,bins=np.linspace(gy.min()*2,gy.max()*2,100))
    b = np.array([np.mean([b[i],b[i+1]]) for i in range(len(b)-1)])
    n2 = gaussian_filter1d(n,1.5/60**2*np.pi/180/2.35/(b[1]-b[0]))
    pdb.set_trace()

    return

def examineSecondBatch(filename):
    """Determine optimal RMS slope in dispersion direction for
    a given measurement
    """
    #Load image
    i,phase,lat = zygo.readzygo(filename)
    lat = lat*1000.
    #Fill missing values with interpolated value
    phase = man.nearestNaN(phase)

    #Loop through rotation angles and compute
    #RMS slope over region of interest
    fom = []
    ptov = []
    for ang in np.linspace(0.,360.,180):
        #Rotate image by appropriate angle
        d = rotate(phase,ang)
        #Select region of interest
        sh = np.shape(d)
        d = d[round(sh[0]/2-25./2/lat):round(sh[0]/2+25./2/lat),\
                round(sh[1]/2+(47.-32./2)/lat):round(sh[1]/2+(47.+32./2)/lat)]
        #Remove tilts
        res = fit.legendre2d(d,xo=1,yo=1)
        d = d-res[0]
        #Get ptov fom
        ptov.append(anal.ptov(d))
        #Transform to slope in dispersion direction
        sl = np.diff(d*1000./lat*180/np.pi*60**2,axis=0)
        sl2 = np.copy(sl)
##        sl[round(sh[0]/2-25./2/lat):round(sh[0]/2+25./2/lat),\
##                round(sh[1]/2+(42.-32./2)/lat):round(sh[1]/2+(42.+32./2)/lat)]\
##                = np.nan

        
##        sl = sl[round(sh[0]/2-25./2/lat):round(sh[0]/2+25./2/lat),\
##                round(sh[1]/2+(47.-32./2)/lat):round(sh[1]/2+(47.+32./2)/lat)]

        #Compute rms of slope
        sl = np.abs(sl - np.mean(sl)) #Subtract average tilt
        sl = sl.flatten()
        sl = np.sort(sl)
        #Histogram slopes
##        y,x = np.histogram(sl,bins=np.linspace(-200.,200.,1000))
##        pdb.set_trace()
        #fom.append(anal.fwhm(x[1:],y))
        fom.append(sl[round(.875*np.size(sl))]-sl[round(.125*np.size(sl))])
    return np.array(fom),np.array(ptov)
