import numpy as np
import matplotlib.pyplot as plt
import pyfits
import utilities.imaging.stitch as st
import zernikemod as zern
import pdb

def analyzeRepeatability(file1,file2):
    """
    Read in metrology data, overlap the images, fit to
    first 20 Zernike terms, return repeatability image
    """
    #Load data
    img1 = pyfits.getdata(file1)
    img2 = pyfits.getdata(file2)

    #Overlap images
    img2 = st.overlapImages(img1,img2)
    resid = img1-img2

    #Perform Zernike fit
    cx,cy,rad = zern.locateimage(resid,35.,45.)
    coeff = zern.zcoeff(resid,cx=cx,cy=cy,rad=rad)

    #Return filtered repeatability
    return coeff[-1]

def analyzeFirstGratings(fname):
    """Read in the figure data from WFS measurements of
    the trilayer coated gratings. Select the region of interest
    assuming the grating is in the center, and add 2 mm of
    margin to perimeter.
    (38 x 30)
    Compute gradient in dispersion direction and compute
    76% encircled slope figure of merit.
    """
    #Load data
    img = pyfits.getdata(fname)

    #Select region of interest
    img = img[64-39/2:64+39/2,64-16:64+16]

    #Compute gradient
    dx = 76.2/86.25 #mm per pix
    gy,gx = np.gradient(img,dx*1000)

    #Bin up and get fom
    pdb.set_trace()
    gx = gx - np.nanmean(gx)
    gx = gx.flatten()*180/np.pi*60**2
    gx = gx[~np.isnan(gx)]
    gx.sort()
    pdb.set_trace()
    return gx[round(.875*np.size(gx))]-gx[round(.125*np.size(gx))]
