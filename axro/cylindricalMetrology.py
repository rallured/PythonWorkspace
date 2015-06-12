import numpy as np
from midfreq import removecyl
import glob
import matplotlib.pyplot as plt
import legendremod as leg
import pdb
import scipy.optimize

def constructCyl(piston,tip,tilt,radius,ang,sh=(128,128)):
    """Construct a 2D cylinder array with arbitrary radius
    and angle."""
    x,y = np.meshgrid(np.arange(0.,sh[1]),np.arange(0.,sh[0]))
    x = x - np.nanmean(x)
    y = y - np.nanmean(y)
    x2 = x*np.cos(ang) - y*np.sin(ang)
    y2 = x*np.sin(ang) + y*np.cos(ang)
    cyl = (radius/np.abs(radius))*np.sqrt(radius**2-x2**2)
    cyl = cyl - np.nanmean(cyl)
    piston = np.zeros(sh) + piston
    tip = x/x.max()*tip
    tilt = y/y.max()*tilt
    
    return (cyl + piston + tip + tilt)

def residFun(img,sh,p):
    return np.nansum((img-constructCyl(*p,sh=sh))**2)

def removeAlignment(img,cylaxis=0):
    """Attempting to remove rigorous cylinder from metrology data
    as opposed to the first few 2D Legendres, which leave residual
    higher order deformations"""
    #Ensure first axis is cylindrical
    if cylaxis!=0:
        img = np.transpose(img)

    #Define merit function
    sh = np.shape(img)
##    fun = lambda p: np.nansum((img-\
##                            constructCyl(p[0],p[1],p[2],p[3],[4],sh=sh))**2)
    fun = lambda p: residFun(img,sh,p)

    #Estimate initial guesses
    res = leg.leg2dfit(img,2,3)[0]
    start = np.array([np.nanmean(img),res[0,1],res[1,0],\
                      -res[0,2]/2.,1.*np.pi/180])
    #Perform optimization
    res = scipy.optimize.minimize(fun,start,method='nelder-mead',\
                  options={'ftol':1.e-2,'disp':True})

    #Construct fitted alignment term
    align = constructCyl(*res['x'],sh=sh)

    pdb.set_trace()

    return align

def plotClip(clip,fig):
    """Plot the three small aperture measurements for a given clip
    clip is a string 'TopRight' or 'TopMidRight' as per Ben's
    naming convention.
    """
    left = np.genfromtxt('150529'+clip+'_Left.txt')
    center = np.genfromtxt('150529'+clip+'_Center.txt')
    right = np.genfromtxt('150529'+clip+'_Right.txt')
    #Remove alignment terms, up to 1st order legendre
    fit = leg.leg2dfit(left,2,5,reconstruct=True)
    left = left-fit[2]
    fit = leg.leg2dfit(center,2,5,reconstruct=True)
    center = center-fit[2]
    fit = leg.leg2dfit(right,2,5,reconstruct=True)
    right = right-fit[2]

##    fig = plt.figure()
    fig.clf()
    fig.add_subplot(131)
    plt.imshow(np.fliplr(left))
    plt.colorbar()
    fig.add_subplot(132)
    plt.imshow(np.fliplr(center))
    plt.colorbar()
    plt.title(clip)
    fig.add_subplot(133)
    plt.imshow(np.fliplr(right))
    plt.colorbar()
