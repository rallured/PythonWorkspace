import numpy as np
from matplotlib.pyplot import *
from os import *

def createarrays(scanfile):
    #Read data in from file, columns are xpos, ypos, and centroid
    data = np.genfromtxt(scanfile,skip_header=1)
    data = np.transpose(data)
    xpos = data[0]
    ypos = data[1]
    cent = data[2]

    #Create 2D arrays for the contour plot
    xuniq = np.unique(data[0])
    yuniq = np.unique(data[1])
    xarr = np.zeros((size(xuniq),size(yuniq)))
    yarr = np.zeros((size(xuniq),size(yuniq)))
    carr = np.zeros((size(xuniq),size(yuniq)))

    #Populate arrays
    for i in range(size(xpos)):
        xind = np.where(xuniq == xpos[i])
        yind = np.where(yuniq == ypos[i])
        xarr[xind,yind] = xpos[i]
        yarr[xind,yind] = ypos[i]
        carr[xind,yind] = cent[i]

    return (xarr,yarr,carr)

def scancontour(xarr,yarr,carr):
    #Set up plot window
    clf()
    hold(True)
    
    lev = np.arange(0,np.max(carr)+20,int(np.max(carr))/40.)
    print lev
    contourf(xarr,yarr,carr,levels=lev)
    colorbar()
    title('Uniformity Scan')
