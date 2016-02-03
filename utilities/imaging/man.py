#This submodule includes routines to manipulate image arrays
from numpy import *
from scipy.interpolate import griddata

import pdb

def unpackimage(data,xlim=[-1,1],ylim=[-1,1],remove=True):
    """Convert a 2D image into x,y,z coordinates.
    x will relate to 2nd index in order to correspond to abscissa in imshow
    y will relate to 1st index in order to correspond to oordinate in imshow
    if remove is True, NaNs will not be returned in the list of coordinates
    """
    x,y = meshgrid(linspace(xlim[0],xlim[1],shape(data)[1]),\
                   linspace(ylim[0],ylim[1],shape(data)[0]))
##    y = y.flatten()
##    x = x.flatten()
    
    if remove is True:
        ind = invert(isnan(data.flatten()))
        return x.flatten()[ind],y.flatten()[ind],data.flatten()[ind]

    return x.flatten(),y.flatten(),data.flatten()

def shiftNaN(img,n=1,axis=0):
    """This function shifts an image in a NaN padded array
    Specify which axis to shift, and specify which direction
    """
    #Construct array to insert
    if axis is 0:
        ins = repeat(nan,abs(n)*shape(img)[1]).reshape(abs(n),shape(img)[1])
    else:
        ins = repeat(nan,abs(n)*shape(img)[0]).reshape(abs(n),shape(img)[0])
    #If direction=0, shift to positive
    if n > 0:
        img = delete(img,arange(shape(img)[1]-n,shape(img)[1]),axis=axis)
        img = insert(img,0,ins,axis=axis)
    else:
        n = abs(n)
        img = delete(img,arange(n),axis=axis)
        img = insert(img,-1,ins,axis=axis)
    return img

def padNaN(img,n=1,axis=0):
    """Pads an image with rows or columns of NaNs
    If n is positive, they are appended to the end of
    the specified axis. If n is negative, they are
    appended to the beginning
    """
    #Construct array to insert
    if axis is 0:
        ins = repeat(nan,abs(n)*shape(img)[1]).reshape(abs(n),shape(img)[1])
    else:
        ins = repeat(nan,abs(n)*shape(img)[0]).reshape(abs(n),shape(img)[0])
        ins = transpose(ins)
    #If direction=0, shift to positive
    if n < 0:
        img = concatenate((ins,img),axis=axis)
    else:
        img = concatenate((img,ins),axis=axis)
    return img

def padRect(img):
    """Pads an image with an outer NaN rectangle"""
    img = padNaN(img,n=1,axis=0)
    img = padNaN(img,n=-1,axis=0)
    img = padNaN(img,n=1,axis=1)
    img = padNaN(img,n=-1,axis=1)
    return img

def tipTiltPiston(img,piston,tip,tilt,tx=None,ty=None):
    """This function adds a constant and
    tip and tilt to an array
    This makes use of tilt arrays tx,ty
    If not provided, compute using meshgrid
    Updated
    """
    if tx is None:
        ty,tx = meshgrid(arange(shape(img)[1]),\
                                arange(shape(img)[0]))
        tx = (tx-mean(tx)) / tx.max()
        ty = (ty-mean(ty)) / ty.max()

    return img + piston + tip*tx + tilt*ty

def nearestNaN(arr,method='nearest'):
    """Fill the NaNs in a 2D image array with the griddata
    nearest neighbor interpolation"""
    ishape = shape(arr)
    #Unpack image both with and without removing NaNs
    x0,y0,z0 = unpackimage(arr,remove=False)
    x1,y1,z1 = unpackimage(arr,remove=True)

    #Interpolate onto x0,y0 grid
    newarr = griddata((x1,y1),z1,(x0,y0),method=method)

    return newarr.reshape(ishape)

def rebin(a,shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return nanmean(nanmean(a.reshape(sh),axis=3),axis=1)

def stripnans(d):
    if len(shape(d)) is 1:
        return d[~isnan(d)]
    newsize = shape(d)[1]
    while sum(isnan(d[0]))==newsize:
        d = d[1:]
    while sum(isnan(d[-1]))==newsize:
        d = d[:-1]
    newsize = shape(d)[0]
    while sum(isnan(d[:,0]))==newsize:
        d = d[:,1:]
    while sum(isnan(d[:,-1]))==newsize:
        d = d[:,:-1]
    return d

def transformation(x,y,r=0.,tx=0.,ty=0.):
    """Return x and y vectors after applying a rotation about
    the origin and then translations in x and y
    """
    x,y = cos(r)*x+sin(r)*y,-sin(r)*x+cos(r)*y
    x,y = x+tx,y+ty
    return x,y
