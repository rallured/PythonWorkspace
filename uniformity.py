import os
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import gaussfitter as gf
import curvefitting

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
    x=xy[:,0]
    y=xy[:,1]
    plt.clf()
    plt.plot(y)
    plt.xlim(0, 512)
    plt.show()
    if raw_input('Discard Spectrum (Y,y) ? ') in ('Y', 'y'):
        centroid=np.zeros(2)
	chisquared = 0
	print centroid
    else:
   	#x0=raw_input('Enter guess for initial x ')
    	#x0=float(x0)
    	#xm=raw_input('Enter guess for max x ')
    	#xm=float(xm)
    	#print x0, xm
	plt.clf()
 	plt.plot(xy[:,1])
 	plt.xlim(0,512)
	plt.show()
	pars=curvefitting.gaussguess(x[80.:200.], y[80.:200.],nterms=4)
	print ('fit parameters  = '), pars
	params = gf.onedgaussfit(x[80.:200.], y[80.:200.], params=pars)
	plt.plot(x[80.:200.],params[1])
   
    	while (raw_input('ReFit (y)? ') == 'y'):

		x0=raw_input('Enter new initial x ')
		x0=float(x0)
		xm=raw_input('Enter new max x ')
		xm=float(xm)
		plt.clf()
 		plt.plot(xy[:,1])
 		plt.xlim(0,512)
		plt.show()
		pars=curvefitting.gaussguess(x[x0:xm], y[x0:xm],nterms=4)
		print ('fit parameters  = '), pars
		params = gf.onedgaussfit(x[x0:xm], y[x0:xm], params=pars)
		plt.plot(x[x0:xm],params[1])

    	centroid = params[0][2]
    	chisquared = params[3]
    
    return (centroid, chisquared, xy[:,0], xy[:,1])

    
def parse(dir = './'):
    centroids = []
    
    filelist = os.listdir(dir)
    
    #These next four lines assume a filename of the form:
    #'specx_[0-9][0-9]y_[0-9][0-9]'
    xmin = int(filelist[0].split('_')[1][:-1])
    ymin = int(filelist[0].split('_')[2])
    xmax = int(filelist[len(filelist)-1].split('_')[1][:-1])
    ymax = int(filelist[len(filelist)-1].split('_')[2])
    
    for fname in filelist:
##        x = fname.split('_')[1][:-1]
##        y = fname.split('_')[2]
        
        print('Currently on spectrum: ' + fname)
        centroid_candidate, chisq, x, y = auto_fit(fname)
        print centroid_candidate
        #Give the user the option to mark this image bad and move on.
        print chisq
	if chisq == 0:
		centroids.append(0)
	elif raw_input('Discard Fit? ') in ('Y', 'y'):
	    #Ask for bounds and replot, bounds ok?

	    #Ask for initial fit parameters, use these to call gaussfit

	    #Fit ok?
            	centroids.append(0)

        else:
		centroids.append(centroid_candidate)
    
    xsteps = np.arange(xmin, xmax+1)
    ysteps = np.arange(ymin, ymax+1)
    zvals  = np.array(centroids).reshape(xmax-xmin+1, ymax-ymin+1).transpose()
    
    savefile = open('FieldUniformityData.dat', 'wb')
    np.savetxt(savefile, zvals, fmt='%4.7f', delimiter='\t')#, newline='\n')
    savefile.close()

    plt.figure(2)
    plt.clf()
    plt.contourf(xsteps, ysteps, zvals, np.arange(75, 101, 3))
    plt.title('Field Uniformity Scan')
    plt.colorbar()
    plt.savefig("FieldUniformScan.png")
    plt.show()
    
    return

    
def replot(filename = 'FieldUniformityData.dat', xmin=21, ymin=7, xmax=39, ymax=23):
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

def addpositionkeywords(dir='./'):
    filelist = os.listdir(dir)
    print filelist

    for fname in filelist:
        #Get x position and y position
        spl = fname.split('_')
        x = spl[1][0:2]
        x = int(x)
        y = spl[2]
        y = int(y)
        
        #Open fits file and add keywords
        f = pyfits.open(fname,mode="update")
        head = f[1].header
        head.update('XPOS',x)
        head.update('YPOS',y)

        #Saves the altered fits file
        f.close()

