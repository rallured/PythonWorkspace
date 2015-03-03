from numpy import *
from matplotlib.pyplot import *
import legendremod as leg
import pdb
import reconstruct as rec
from plotting import nanmean

#Get power distribution array from Thorlabs CSV
def getPow(filename):
    if filename.split('.')[1]=='txt':
        return genfromtxt(filename)
    
    f = open(filename,'r')

    l = f.readline()
    while l != '*** POWER DISTRIBUTION [a.u.] ***\r\n':
        l = f.readline()
    l = f.readline()
    px = int(l.split(',')[1])
    l = f.readline()
    py = int(l.split(',')[1])

    power = zeros((px,py))
    for i in range(py):
        l = f.readline()
        l = l.split(',')[:-1]
        l = array(l)
        l = l.astype('float')
        power[:,i] = l

    return power

#Get spot positions from Thorlabs CSV
def getSpots(filename):
    f = open(filename,'r')

    l = f.readline()
    while l != '*** CENTROIDS [pixels] ***\r\n':
        l = f.readline()
    l = f.readline()
    px = int(l.split(',')[1])
    l = f.readline()
    py = int(l.split(',')[1])

    cent = zeros((px*2,py))
    for i in range(py):
        l = f.readline()
        l = l.split(',')[:-1]
        l = array(l)
        l = l.astype('float')
        cent[:,i] = l

    x = cent[0::2,:]
    y = cent[1::2,:]

    return x,y

#Get wavefront measurement from power distribution, spot positions
#and reference spot positions
def getWave(figspots,figpow,refspots):
    #Get spots
    x,y = getSpots(figspots)
    xr,yr= getSpots(refspots)
    #Get power
    p = getPow(figpow)

    #Set up NaNs
    ind = isnan(x)
    xr[ind] = NaN
    yr[ind] = NaN
    p[ind] = NaN
    ind = isnan(xr)
    x[ind] = NaN
    y[ind] = NaN
    p[ind] = NaN
    ind = where(p==0.)
    x[ind] = NaN
    y[ind] = NaN
    xr[ind] = NaN
    yr[ind] = NaN

    #Set up arrays for reconstruction
    ind = isnan(x)
    xs = (x-xr)/5200.
    ys = (y-yr)/5200.
    xs = xs - nanmean(xs)
    ys = ys - nanmean(ys)
    xs[ind] = 100.
    ys[ind] = 100.
    ph = zeros(shape(xs))
    ph[ind] = 100.
    pdb.set_trace()

    #Reorder arrays as Fortran
    xs = asarray(xs,order='F')
    ys = asarray(ys,order='F')
    ph = asarray(ph,order='F')

    #Reconstruct and return wavefront
    wv = rec.reconstruct(xs,ys,1.e-12,150.,ph)
    wv[wv==100.] = NaN

    return wv
