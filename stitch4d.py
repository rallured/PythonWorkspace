from numpy import *
from matplotlib.pyplot import *
from scipy.interpolate import griddata
from legendremod import unpackimage
import pdb
from plotting import nanmean
import scipy.optimize

###This function will add translation, piston, and tilt to array2,
###interpolate to array1, and calculate rms difference
##def stitchArrays(array1,array2,tx,ty,p,ax,ay):
##    #Add piston, tip, and tilt to array2
##    x,y = meshgrid(arange(982.),arange(1002.))
##    array2 = array2 + p + x*ax + y*ay
##    
##    #Unpack arrays into x,y,z lists
##    x2,y2,z2 = unpackimage(array2)
##    x1,y1,z1 = unpackimage(array1)
##
##    #Add translations
##    x2 = x2 + tx
##    y2 = y2 + ty
##
##    print 'interpolate'
##    #Interpolate onto array1
##    iz = griddata((x2,y2),z2,(x1,y1))
##
##    pdb.set_trace()
##
##    return

x,y = meshgrid(arange(982.),arange(1002.))

#This function shifts an image in a NaN padded array
#Specify which axis to shift, and specify which direction
def shiftNaN(img,n=1,axis=0):
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

#Transform array
def transformArray(array,ty,tx,p,ax,ay):
    array = shiftNaN(array,n=round(ty),axis=0)
    array = shiftNaN(array,n=round(tx),axis=1)
    array = array + p + x*ax + y*ay
    return array

#Instead assume integer y translations
def overlapArrays(array1,array2,ty,tx,p,ax,ay):
    array2 = transformArray(array2,ty,tx,p,ax,ay)

    #Compute RMS diff
    rms = sqrt(nanmean((array2-array1)**2))

    print ty,tx,p,ax,ay,rms
        
    return rms

def matchArrays(array1,array2,ty,tx,pis,ax,ay):
    #Create function
    fun = lambda p: overlapArrays(array1,array2,p[0],p[1],p[2],p[3],p[4])

    #Optimize function
    start = zeros(5)
    start[0] = ty
    start[1] = tx
    start[2] = pis
    start[3] = ax
    start[4] = ay
    res = scipy.optimize.minimize(fun,start,method='nelder-mead',\
                                  options={'disp':True,'maxfev':1000})
    print res['x']

    return transformArray(array2,*res['x'][:])

