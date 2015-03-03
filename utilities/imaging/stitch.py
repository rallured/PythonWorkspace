#This submodule facilitates stitching of metrology images
#The images must have fiducials to compute translation and rotation
#Piston/tip/tilt is done by arrays after translation/rotation is fixed
#by fiducials
import utilities.transformations as tr
from numpy import *
from matplotlib.pyplot import *
import man
from scipy.optimize import minimize
from analysis import getPoints
import pdb
from scipy.interpolate import griddata
from utilities.plotting import nanmean

def transformCoords(x,y,tx,ty,theta):
    """Transforms coordinates x,y by translating tx,ty
    and rotation theta about x
    Returns: x,y of new coords
    """
    trans = tr.translation_matrix([tx,ty,0])
    rot = tr.rotation_matrix(theta,[0,0,1],point=[mean(x),mean(y),0])
    pos0 = array((x,y,repeat(0.,size(x)),repeat(1.,size(x))))
    pos1 = dot(trans,dot(rot,pos0))

    return pos1[0],pos1[1]

def sumOfSquares(x1,y1,x2,y2):
    """Computes the sum of the squares of the residuals
    for two lists of coordinates
    """
    return sum(sqrt((x1-x2)**2+(y1-y2)**2))
    

def matchFiducials(x1,y1,x2,y2):
    """This function will compute a rotation and translation
    to match a list of fiducial coordinates
    Returns: translation tx,ty and rotation theta about zhat
    to bring x2,y2 to x1,y1
    """
    #Make function to minimize
    fun = lambda p: sumOfSquares(x1,y1,*transformCoords(x2,y2,*p))

    #Make starting guess
    start = zeros(3)
    start[0] = mean(x1-x2)
    start[1] = mean(y1-y2)
    start[2] = .0001

    #Run minimization and return fiducial transformation
    res = minimize(fun,start,method='nelder-mead',\
                   options={'disp':True,'maxfev':1000})
    
    return res['x']

def matchPistonTipTilt(img1,img2):
    """This function applies piston and tip/tilt
    to minimize RMS difference between two arrays
    Returns: img2 matched to img1
    """
    #Make function to minimize
    fun = lambda p: nanmean((img1 - man.tipTiltPiston(img2,*p))**2)

    #Run minimization and return matched image
    res = minimize(fun,[1.,.1,.1],method='nelder-mead',\
                   options={'disp':True,'maxfev':1000})

    return man.tipTiltPiston(img2,*res['x'])

def stitchImages(img1,img2):
    """Allows user to pick fiducials for both images.
    Function then computes the transform to move img2
    to img1 reference frame.
    Updated
    """
    #Get both fiducials
    xf1,yf1 = getPoints(img1)
    xf2,yf2 = getPoints(img2)

    #Match them
    tx,ty,theta = matchFiducials(xf1,yf1,xf2,yf2)

    #Pad img1 based on translations
    img1 = man.padNaN(img1,n=round(tx),axis=1)
    img1 = man.padNaN(img1,n=round(ty),axis=0)
    #Shift img1 fiducial locations
    if tx<0:
        xf1 = xf1 - tx
    if ty<0:
        yf1 = yf1 - ty

    #Get x,y,z points from stitched image
    x2,y2,z2 = man.unpackimage(img2,xlim=[0,shape(img2)[1]],\
                           ylim=[0,shape(img2)[0]])

    #Apply transformations to x,y coords
    x2,y2 = transformCoords(x2,y2,ty,tx,theta)

    #Get x,y,z points from reference image
    x1,y1,z1 = man.unpackimage(img1,remove=False,xlim=[0,shape(img1)[1]],\
                           ylim=[0,shape(img1)[0]])

    #Interpolate stitched image onto expanded image grid
    newimg = griddata((x2,y2),z2,(x1,y1),method='linear')
    print 'Interpolation ok'
    newimg = newimg.reshape(shape(img1))

    #Images should now be in the same reference frame
    #Time to apply tip/tilt/piston to minimize RMS
    newimg = matchPistonTipTilt(img1,newimg)

    #Would like list of enlarge image showing all valid data, this is final step
    #Avoid overwritting fiducials
    #Save indices of NaNs near fiducials
    find = logical_and(sqrt((y1-xf1[0])**2+(x1-yf1[0])**2) < 15.,\
                       isnan(img1).flatten())
    for i in range(1,size(xf1)):
        find = logical_or(find,\
                logical_and(sqrt((y1-xf1[i])**2+(x1-yf1[i])**2) < 15.,\
                        isnan(img1).flatten()))

    #Include newimg data
    img1[isnan(img1)] = newimg[isnan(img1)]
    #Reset fiducials to NaNs
    img1[find.reshape(shape(img1))] = NaN

    #Return translations to automatically pad next image
    return img1,tx,ty
