import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import utilities.imaging.analysis as anal
import scipy.interpolate as interp
from scipy.ndimage import rotate
import pdb,sys

def dispersion(order):
    """Return the linear dispersion of the test for a given order
    """
    L = 8052.e-3 #Distance from grating to detector plane
    d = 160.e-9 #Period of grating grooves
    i = 1.5*np.pi/180 #Incidence angle of grating
    return order*L/d/np.cos(i)

def nominalX(energy,order):
    """Return the nominal x displacement for a given order
    of diffraction for a given wavelength.
    Give energy in eV
    Returns x in mm
    """
    wave = 1239.8 / energy
    return dispersion(order)*(wave*1e-6)

def linesep(energy1,order1,energy2,order2):
    """Return the distance between two lines in mm using the
    theoretical dispersion."""
    return (nominalX(energy1,order1)-nominalX(energy2,order2))

def computeHEW(h):
    """Compute the HEW assuming negligible background.
    Image has already been collapsed in proper axis.
    Simply compute CDF, interpolate with a spline, and
    return the HEW
    """
    cdf = np.cumsum(h)/np.sum(h) #Normalized CDF
    pix = np.arange(np.size(cdf)) #X vector
    sp = interp.interp1d(pix,cdf,kind=3) #Cubic spline fit
    pix = np.linspace(pix[0],pix[-1],np.size(cdf)*10.) #Resampled x vector
    cdf = sp(pix) #Resampled CDF
    return pix[np.argmin(np.abs(cdf-.75))] - \
           pix[np.argmin(np.abs(cdf-.25))]

#Processing steps
##1. Select subaperture - by clicking
##2. Apply rotation, if any
##3. Collapse in cross-dispersion direction
##4. Cumsum and spline to create CDF
##5. Return HEW
##6. Compute uncertainty a la DeRoo, need to convert to electrons first
def analyzePANTER(img):
    #Clean hot pixels
    fig = plt.figure()
    plt.imshow(img,norm=LogNorm())
    plt.colorbar()
    flag = False
    while flag is False:
        mx = input('Max pixel value?')
        img[img>mx] = 0.
        plt.clf()
        plt.imshow(img,norm=LogNorm())
        plt.colorbar()
        flag = input('OK? (0/1)')
    plt.close(fig)
    
    #Get subaperture by clicking two points on the image
    x,y = anal.getPoints(img,log=True)
    img = img[y.min():y.max(),\
              x.min():x.max()]

    #Loop through rotation angles and tabulate HEW vs. angle
    hew = np.zeros(101)
    r = np.linspace(-2.,2.,101)
    for i in range(np.size(r)):
        newimg = rotate(img,r[i],order=1) #Linear spline rotation
        sys.stdout.write('Rotation Done\r')
        hew[i] = computeHEW(np.sum(newimg,axis=0)) #Compute HEW
        sys.stdout.write('HEW Done %03i\r' % i)

    #Plot the CDF
    newimg = rotate(img,r[np.argmin(hew)],order=1)
    fig = plt.figure()
    plt.plot(np.sum(newimg,axis=0))
    fig2 = plt.figure()
    plt.plot(np.sum(newimg,axis=1))

    #Return the list of HEWs
    return hew,newimg
