import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import utilities.imaging.analysis as anal
import scipy.interpolate as interp
import scipy.ndimage as nd
import pdb,sys,pyfits,glob
import scipy.signal

arcusdir = '/home/rallured/Dropbox/Arcus/Alignment/PANTER/data/'

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
    h = h.astype('float')
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
def analyzePANTER(img,minang=-2.,maxang=2.):
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
    r = np.linspace(minang,maxang,101)
    for i in range(np.size(r)):
        newimg = nd.rotate(img,r[i],order=1) #Linear spline rotation
        sys.stdout.write('Rotation Done\r')
        hew[i] = computeHEW(np.sum(newimg,axis=0)) #Compute HEW
        sys.stdout.write('HEW Done %03i\r' % i)

    #Plot the CDF
    newimg = nd.rotate(img,r[np.argmin(hew)],order=1)
    fig = plt.figure()
    plt.plot(np.sum(newimg,axis=0))
    fig2 = plt.figure()
    plt.plot(np.sum(newimg,axis=1))

    #Return the list of HEWs
    return hew,newimg

def estimateYawAlign():
    #Load in Tropic data
    d60 = pyfits.getdata(arcusdir+'10272014/Tropic/HK141027.060_img.fits')
    d74 = pyfits.getdata(arcusdir+'10272014/Tropic/HK141027.074_img.fits')

    #Get subapertured images collapsed
    ref = np.sum(d60[:400,350:500],axis=1)
    act = np.sum(d60[400:,350:500],axis=1)
    ref = ref.astype('float')
    act = act.astype('float')
    
    #Get initial centroids of reference vs. active
    actf = scipy.signal.savgol_filter(act,101,8,mode='constant')
    maxpix = actf.argmax()
    centroid = np.sum(act/np.sum(act)*np.arange(560.))
    centoffset = centroid - maxpix
    refcentroid = np.sum(ref/np.sum(ref)*np.arange(400.))
    
    #Find difference of peak to centroid in active data
    final = np.sum(d74,axis=1)
    finalf = scipy.signal.savgol_filter(final,101,8,mode='constant')
    peak = finalf.argmax()
    actcent = peak + centoffset

    return actcent,refcentroid

def pitchScan():
    #Load in Tropic data
    fn = glob.glob(arcusdir+'10272014/Tropic/HK*.fits')
    fn.sort()
    fn = fn[5:14] #Select pitch scan

    #Loop through and compute RMS figure of merit
    fom = []
    for f in fn:
        d = pyfits.getdata(f)
        d = np.sum(d,axis=1)
        x = np.arange(np.size(d))
        cent = np.average(x,weights=d)
        fom.append(np.sqrt(np.average((x-cent)**2,weights=d)))

    return fom

def rollScan():
    #Load in Tropic data
    fn = glob.glob(arcusdir+'10272014/Tropic/HK*.fits')
    fn.sort()
    fn = fn[14:22]
    order = [0,1,2,3,4,7,6,5]
    fn = [fn[i] for i in order]

    #Loop through and compute RMS figure of merit
    fom = []
    for f in fn:
        d = pyfits.getdata(f)
        d = np.sum(d[:,400:600],axis=0)
        x = np.arange(np.size(d))
        cent = np.average(x,weights=d)
        fom.append(np.sqrt(np.average((x-cent)**2,weights=d)))

    return fom

def zeroOrder():
    #Load in Tropic zero order file
    d = pyfits.getdata(arcusdir+'10242014/Tropic/HK141027.050_img.fits')

def mosaic():
    #Load in PSPC mosaic file
    mosaic = pyfits.getdata(arcusdir+'10242014/ARCUS_Grating_Mosaic_MgK.fits')

    fig = plt.figure()
    plt.imshow(np.flipud(np.transpose(mosaic)),norm=LogNorm())
    plt.title('PSPC Mosaic Image')
    plt.plot([1500.,1500.],[0.,1483.],'k--')
    plt.ylim([1483.,0.])

    fig = plt.figure()
    plt.semilogy(np.sum(mosaic[1500:],axis=1))
    plt.text(120.,.04,'Mg-K 1')
    plt.text(310.,.1,'Cu-L 1')
    plt.text(410.,.04,'Co-L 1')
    plt.text(600.,.3,'Fe-L 1')
    plt.text(900.,4.5,'O-K 1')
    plt.text(1075.,.1,'Cu-L 2')
    plt.title('Collapsed PSPC Mosaic Spectrum')
    plt.xlabel('PSPC Pixel')
    plt.ylabel('Counts')
    
def Mg1Lines():
    #Load in PIXI data
    d41 = pyfits.getdata(arcusdir+'10272014/PIXI/IM141027_41.fits')
    d2 = pyfits.getdata(arcusdir+'10292014/IM141029_2.fits')
    
