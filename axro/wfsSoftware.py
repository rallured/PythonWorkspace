import numpy as np
import matplotlib.pyplot as plt
import utilities.imaging.analysis as anal
import os,glob,pyfits

def getSlopes():
    """
    Grab SlopesX and SlopesY in current directory
    """
    filenames = glob.glob('SlopesX*')
    filenames.sort()
    slopesx = [pyfits.getdata(s) for s in filenames]
    rmsx = [anal.rms(s) for s in slopesx]
    filenames = glob.glob('SlopesY*')
    filenames.sort()
    slopesy = [pyfits.getdata(s) for s in filenames]
    rmsy = [anal.rms(s) for s in slopesy]

    return slopesx,rmsx,slopesy,rmsy

if __name__=='__main__':
    #Investigate repeatability data
    os.chdir('/home/rallured/Dropbox/WFS/SystemAlignment/Repeatability/SDK/Single')
    sx_s1,rx_s1,sy_s1,ry_s1 = getSlopes()
    os.chdir('../100')
    sx_s100,rx_s100,sy_s100,ry_s100 = getSlopes()
    os.chdir('../../GUI/Single')
    sx_g1,rx_g1,sy_g1,ry_g1 = getSlopes()
    os.chdir('../100')
    sx_g100,rx_g100,sy_g100,ry_g100 = getSlopes()
