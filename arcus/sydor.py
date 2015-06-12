import numpy as np
from scipy.ndimage import rotate
from zygo import readzygo
import pdb,glob
import matplotlib.pyplot as plt

#Need to analyze Sydor metrology of glass wafers
#Load in metrology data, and compute slope distributions
#over 95 x 75 mm rectangle about the center of the
#image. Scan through a rotation angle. NaNs filled in
#with mean of image. Should be small effect on slope
#statistics. Image shape will be changed by rotation, so
#choose center of image independently of shape.

def slopeStatistics(filename):
    #Load image
    i,phase,lat = readzygo(filename)
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
        fom.append(sl[round(.875*np.size(sl))]-sl[round(.125*np.size(sl))])
    return fom

#Loop through Sydor ascii files and plot FoM
#for each one
def sydorPlot():
    #Get filenames
    f = glob.glob('/home/rallured/Dropbox/Arcus/Sydor/*.asc.txt')
    #Loop through and plot FoM vs. rotation
    ang = np.linspace(0.,180.,180)
    plt.ion()
    plt.figure()
    for fi in f:
        plt.plot(ang,slopeStatistics(fi),label=fi.split('.')[0])
    #Plot requirement
    plt.([0.,180],[17.7,17.7],'k--')
    plt.title('Gaussian FWHM Equivalent Widths')
    plt.xlabel('Wafer Rotation Angle')
    plt.ylabel('Width (arcsec)')
    
