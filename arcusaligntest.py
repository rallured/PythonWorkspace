from numpy import *
from matplotlib.pyplot import *
import tifffile as tiff
import sys
from scipy.interpolate.interpolate import interp1d

#Read in list of images and average their pixel values
def avgImgList(flist):
    im = PIL.Image.open(flist[0])
    p = array(im)
    paccum = copy(p)
    for f in flist[1:]:
        im = PIL.Image.open(f)
        p = array(im)
        paccum = paccum + p
    return paccum / size(flist)

#Load in sequence, get rid of redundant data
def loadSequence(fname):
    p = tiff.imread(fname)
    p = p[:,:,:,0]
    return p

#After establishing xrange for line, interpolate other image
#onto appropriate x vector from reference image
def interpolateLine(xact,yact,xvec):
    #Interpolation function
    interFn = interp1d(xact,yact,kind='linear')
    #Return interpolated yvalues
    return interFn(xvec)

#Apply offset to line data, interpolate,
#and compute cross-correlation coefficient
def crossCorrelation(xref,yref,xact,yact,off):
    ynew = interpolateLine(xact-off,yact,xref)
    return sum(yref*ynew)/sum(yref**2)

#Compute cross correlation coefficient for an array of offsets
def crossOffsets(xref,yref,xact,yact,off):
    coef = []
    for o in off:
        coef.append(crossCorrelation(xref,yref,xact,yact,o))
    return array(coef)

#Compute cross correlation coefficients and find max offset
def optimalOffset(xref,yref,xact,yact,off):
    c = crossOffsets(xref,yref,xact,yact,off)
    fit = polyfit(off,c,2)
    plot(off,c)
    plot(off,polyval(fit,off))
    return off[argmax(polyval(fit,off))]

#Subtract background of image and compute centroid
def computeCentroid(img):
    bg = min(img[600:,1100])
    x = arange(1024.)
    return sum(x[200:340]*(img[200:340,1100]-bg))/sum(img[200:340,1100]-bg)

#Load all the sequences
def loadAll():
    pos = []
    for i in range(1,7):
        d = loadSequence('Sequence'+str(i)+'.tif')
        d = mean(d,axis=0)
        pos.append(d)
    return pos

#Compute shifts with respect to first position
def computeShifts():
    pos = loadAll()
    xref = arange(0.,1024.)*5.24
    yref = pos[0][:,970]
    ind = logical_and(xref>1000.,xref<1800.)

    shift = []
    for i in range(1,6):
        shift.append(optimalOffset(xref[ind],yref[ind],xref,pos[i][:,970],\
                                   linspace(-50.+25*(i-1),50.+25*(i-1),970)))
    return shift

#Make plot for proposal
def propPlot():
    p = array([165.442,191.332,215.589,241.623,266.041,290.946])
    p = p-p[0]
    p = p[1:]
    psig = array([.059,.102,.087,.143,.147,.072])
    s = computeShifts()
    figure()
    plot(p,s,'*',markersize=10)
    plot(p,p-5,'k--')
    plot(p,p+5.,'k--')
    title('SPO Diffraction Pattern Translation Test')
    xlabel('Keyence Probe Translation Measurement (microns)')
    ylabel('Camera Translation Measurement (microns)')
    
    
    
    
