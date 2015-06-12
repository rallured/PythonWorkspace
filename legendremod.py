#This module is a Legendre polynomial fitter to a rectangular image
#Basic approach is to assign each element in the array an x and y coordinate,
#where the x coordinate refers to the first index
#The array is reshaped into a 1D vector z, and the x and y vectors will match up
#with their respective elements
#For example, z[0] is located at the point (x[0],y[0])

#The x and y vectors then need to be rescaled such that they go from -1 to 1
#This maintains the complete, orthonormal nature of the Legendre polynomials

#Then, we can construct a matrix A where each column corresponds with a 2D
#Legendre polynomial image that matches to the rescaled x and y coordinates
#A coefficient vector C is what needs to be found to minimize the difference
#of the fit image A*C with the actual image Z

#The algorithm that finds the coefficients to minimize this difference is
#found in scipy.linalg.lstsq

from numpy import *
import numpy.polynomial.legendre as leg
import pdb
import scipy.linalg as lin

#Takes in a rectangular array and outputs x,y,z vectors
#x is identified as the FIRST index, a higher index indicates a more positive
#coordinate
def unpackimage(data,remove=True):
    #Shape function used to create x and y vectors
    #Both span vectors span -1 to 1 with equal spacings
    xspan = linspace(-1,1,shape(data)[0])
    yspan = linspace(-1,1,shape(data)[1])

    #Loop through data array and assemble x, y, z vectors
    x = zeros(size(data))
    y = copy(x)
    z = copy(x)
    i = 0
    for xi in range(shape(data)[0]):
        for yi in range(shape(data)[1]):
            z[i] = data[xi,yi]
            x[i] = xspan[xi]
            y[i] = yspan[yi]
            i += 1

    #Remove any NaNs
    if remove==True:
        ind = invert(isnan(z))
        z = z[ind]
        x = x[ind]
        y = y[ind]

    return x,y,z

#Need to construct matrix of images for Legendre coefficients
#Feed x and y coordinate vectors as well as maximum order to fit
#in both x and y axes
def imagematrix(x,y,xorder,yorder):
    #Initialize A
    A = zeros((size(x),xorder*yorder))
    xc = zeros(xorder)
    yc = zeros(yorder)
    i = 0
    for xi in range(xorder):
        xc[xi] = 1 #Set current x order
        for yi in range(yorder):
            yc[yi] = 1 #Set current y order
            A[:,i] = leg.legval(x,xc)*leg.legval(y,yc) #Set image for this order
            yc[yi] = 0 #Reset this y order to 0
            i += 1
        xc[xi] = 0 #Reset this x order to 0

    return A

#Perform least squares fit and put coefficients into x,y matrix
#fits up to xorder-1 in x and yorder-1 in y
#for example, setting xorder=yorder=2 fits piston, tip, tilt, and twist,
#where twist is P_1(x)*P_1(y)
def leg2dfit(data,xorder,yorder,reconstruct=False):
    #Create x,y,z fit vectors
    x,y,z = unpackimage(data,remove=False)
    x2,y2,z2 = unpackimage(data,remove=True)

    #Perform fit
    A = imagematrix(x2,y2,xorder,yorder)
    fit = lin.lstsq(A,z2)
    coeff = reshape(fit[0],(xorder,yorder))

    #Compute RMS difference between fit and data
    rms = sqrt(fit[1]/size(z2))

    #Reconstruct image from coefficients
    #If NaNs are involved, the reconstructed image will not be rectangular,
    #so don't reconstruct in that case
    if reconstruct==True:
        A = imagematrix(x,y,xorder,yorder)
        fitimage = zeros(size(data))
        pdb.set_trace()
        for i in range(size(fit[0])):
            fitimage += A[:,i]*fit[0][i]
        fitimage = reshape(fitimage,shape(data),order='C')
        return coeff,rms,fitimage
    else:
        return coeff,rms
