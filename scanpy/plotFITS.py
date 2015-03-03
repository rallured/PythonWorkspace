import os
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import GaussFit as gf

plt.ion()
plt.hold(True)


def importFITS(fname):
    fitsData = pyfits.open(fname)
    tupleArray = fitsData[1].data
    
    x = y = np.zeros(len(tupleArray))
    for i in range(0, len(tupleArray)):
        x[i], y[i] = tupleArray[i]
    
    xy = np.zeros(2*len(tupleArray))
    xy.resize(len(tupleArray),2)
    xy[:,0] = np.arange(0, 512)
    xy[:,1] = y[:]
    
    return(xy)

    
def auto_fit(fname):
    xy = importFITS(fname)
    
    plt.clf()
    plt.plot(xy[:,1])
    plt.xlim(0, 512)
    plt.show()
    
    params = gf.onedgaussfit(xy[:,0], xy[:,1], params=[0,40,100,10])
    plt.plot(params[1])
    centroid = params[0][2]
    chisquared = params[3]
    
    return (centroid, chisquared)

    
def parseData(dir = './'):
    centroids = []
    
    filelist = os.listdir(dir)
    
    #These next four lines assume a filename of the form:
    #'specx_[0-9][0-9]y_[0-9][0-9]'
    xmin = int(filelist[0].split('_')[1][:-1])
    ymin = int(filelist[0].split('_')[2])
    xmax = int(filelist[len(filelist)-1].split('_')[1][:-1])
    ymax = int(filelist[len(filelist)-1].split('_')[2])
    
    for fname in filelist:
#        x = fname.split('_')[1][:-1]
#        y = fname.split('_')[2]
        
        print('Currently on spectrum: ' + fname)
        centroid_candidate, chisq = auto_fit(fname)
        
        #Give the user the option to mark this image bad and move on.
        print chisq
        if raw_input('Discard Fit? ') in ('Y', 'y'):
            centroids.append(0)
        else:
            centroids.append(centroid_candidate)
    
    xsteps = np.arange(xmin, xmax+1)
    ysteps = np.arange(ymin, ymax+1)
    zvals  = np.array(centroids).reshape(xmax-xmin+1, ymax-ymin+1).transpose()
    
    savefile = open('FieldUniformityData.dat', 'wb')
    np.savetxt(savefile, zvals, fmt='%4.7f', delimiter='\t', newline='\n')
    savefile.close()

    plt.figure(2)
    plt.clf()
    plt.contourf(xsteps, ysteps, zvals, np.arange(75, 101, 3))
    plt.title('Field Uniformity Scan')
    plt.colorbar()
    plt.savefig("FieldUniformScan.png")
    plt.show()
    
    return

    
def replot(filename = 'FieldUniformityData.txt', xmin=21, ymin=7, xmax=39, ymax=23):
    savefile = open(filename, 'rb')
    zvals = np.loadtxt(savefile, delimiter='\t')
    savefile.close()

    xsteps = np.arange(xmin, xmax+1)
    ysteps = np.arange(ymin, ymax+1)
    
    plt.figure(2)
    plt.clf()
    plt.contourf(xsteps, ysteps, zvals, np.arange(75, 101, 3))
    plt.title('Field Uniformity Scan')
    plt.colorbar()
    plt.savefig("FieldUniformScan.png")
    plt.show()
    
    
def fixNaming(dir = './'):
    filelist = os.listdir(dir)
    for fname in filelist:
        if int(fname.split('_')[2]) < 10:
            fixedname = fname[:-1] + '0' + fname[-1:]            
            os.rename(fname, fixedname)
    return
    
    
def frange(start, stop, step):
    vals = []
    while start < stop:
        vals.append(start)
        start += step
    return vals
