#This submodule includes routines to manipulate image arrays
from numpy import *

def unpackimage(data,xlim=[-1,1],ylim=[-1,1],remove=True):
    """Convert a 2D image into x,y,z coordinates.
    x will relate to 2nd index in order to correspond to abscissa in imshow
    y will relate to 1st index in order to correspond to oordinate in imshow
    if remove is True, NaNs will not be returned in the list of coordinates
    """
    y,x = meshgrid(linspace(xlim[0],xlim[1],shape(data)[1]),\
                   linspace(ylim[0],ylim[1],shape(data)[0]))

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
    #If direction=0, shift to positive
    if n < 0:
        img = insert(img,0,ins,axis=axis)
    else:
        img = insert(img,-1,ins,axis=axis)
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



    
