import numpy as np
#from numbapro.cudalib import cublas
import pdb
from scipy.optimize import fmin_slsqp
#from numbapro import guvectorize
#from axro.merit import merit
import pyfits
import gc
#import pyOpt

#Need to write merit function calculator for set of influence
#functions and distortion array
#Then, write a function that uses fmin_sqslp to minimize
#sum of square of residuals

#blas = cublas.Blas()
def ampMeritFunction(voltages,distortion,ifuncs):
    """Simple merit function calculator.
    voltages is 1D array of weights for the influence functions
    distortion is 2D array of distortion map
    ifuncs is 4D array of influence functions
    shade is 2D array shade mask
    Simply compute sum(ifuncs*voltages-distortion)**2)
    """
    #Numpy way
    r = np.dot(ifuncs,voltages)-distortion
    res = np.mean((np.dot(ifuncs,voltages)-distortion)**2)
    return res

def ampMeritFunction2(voltages,**kwargs):
    """Simple merit function calculator.
    voltages is 1D array of weights for the influence functions
    distortion is 2D array of distortion map
    ifuncs is 4D array of influence functions
    shade is 2D array shade mask
    Simply compute sum(ifuncs*voltages-distortion)**2)
    """
    #Numpy way
    distortion = kwargs['inp'][0]
    ifuncs = kwargs['inp'][1]
    res = np.mean((np.dot(ifuncs,voltages)-distortion)**2)
    return res, [], 0

def ampMeritDerivative(voltages,distortion,ifuncs):
    """Compute derivatives with respect to voltages of
    simple RMS()**2 merit function
    """
    res = np.dot(2*(np.dot(ifuncs,voltages)-distortion),ifuncs)/\
           np.size(distortion)
    return res

def ampMeritDerivative2(voltages,f,g,**kwargs):
    """Compute derivatives with respect to voltages of
    simple RMS()**2 merit function
    """
    distortion = kwargs['inp'][0]
    ifuncs = kwargs['inp'][1]
    res = np.dot(2*(np.dot(ifuncs,voltages)-distortion),ifuncs)/\
           np.size(distortion)
    return res.tolist(), [], 0

def rawOptimizer(ifs,dist,bounds=None,smin=0.,smax=5.):
    """Assumes ifs and dist are both in slope or amplitude space.
    No conversion to slope will occur."""
    #Create bounds list
    if bounds is None:
        bounds = []
        for i in range(np.shape(ifs)[0]):
            bounds.append((smin,smax))

    #Get ifs in right format
    ifs = ifs.transpose(1,2,0) #Last index is cell number

    #Reshape ifs and distortion
    sh = np.shape(ifs)
    ifsFlat = ifs.reshape(sh[0]*sh[1],sh[2])
    distFlat = dist.flatten()

    #Call optimizer algoritim
    optv = fmin_slsqp(ampMeritFunction,np.zeros(sh[2]),\
                      bounds=bounds,args=(distFlat,ifsFlat),\
                      iprint=2,fprime=ampMeritDerivative,iter=200,\
                      acc=1.e-10)

    #Reconstruct solution
    sol = np.dot(ifs,optv)

    return sol,optv
    

def slopeOptimizer2(dslopes=None,ifslopes=None,ifuncf=None,\
                        distortionf=None,shadef=None,\
                        dx=100./150,azweight=.015,\
                        smin=0.,smax=5.,bounds=None):
    """Format the arrays to matrices and vectors for
    the merit function. Then call fmin_slsqp to determine
    the optimal voltages to correct the mirror.
    All conversion to slope and any frequency filtering must
    occur prior to optimization call.
    Need to calculate required accuracy parameter to get faster
    convergence. Solution speed is on the order of solve_pzt, however!
    Need more RAM for larger optimization problems
    Order more from Dell ASAP
    In the meantime, make use of Tesla K20 for faster merit functions
    Should be able to get it down to much faster than with 8 core processor
    """
    if dslopes is None:
        #Load in distortion, ifuncs, and shade from filenames
        ifuncs = pyfits.getdata(ifuncf)/(dx*1000)*180./np.pi*60**2
        distortion = pyfits.getdata(distortionf)/(dx*1000)*180./np.pi*60**2
        shade = pyfits.getdata(shadef)
        
        #Reshape ifuncs into 3D matrix
        if ifuncs.ndim == 4:
            ishape = np.shape(ifuncs)
            ifuncs = ifuncs.reshape((ishape[0]*ishape[1],ishape[2],ishape[3]))
        ishape = np.shape(ifuncs)
        ifuncs = ifuncs.transpose(1,2,0)#*1.e3 #Get to units of microns
        axif = np.diff(ifuncs,axis=0) #Axial slopes
        axif = axif[:,:-1,:] #Get rid of last column
        azif = np.diff(ifuncs,axis=1) #Azimuthal slopes
        azif = azif[:-1,:,:] #Get rid of last row
        del ifuncs
        gc.collect()

        #Reshape IFs into 2D matrices
        axif = axif.reshape(((ishape[1]-1)*(ishape[2]-1),ishape[0]))
        azif = azif.reshape(((ishape[1]-1)*(ishape[2]-1),ishape[0]))

        #filter ifuncs array with shade mask
        shade = shade[:-1,:-1] #Get rid of last row and column
        shade = shade.astype('bool')
        shade = shade.flatten()
        axif = axif[shade,:]
        azif = azif[shade,:]*azweight

        #Merge slopes into 2D optimization matrix
        ifslopes = np.concatenate((axif,azif),axis=0)

        #Free up some memory
        del axif
        del azif
        gc.collect()

        #Create distortion slope arrays
        axd = np.diff(distortion,axis=0)
        axd = axd[:,:-1]
        azd = np.diff(distortion,axis=1)
        azd = azd[:-1,:]
        del distortion
        gc.collect()

        #FLatten distortion arrays
        axd = axd.flatten()
        azd = azd.flatten()
        axd = axd[shade]
        azd = azd[shade]*azweight

        #Create 1D distortion slope target
        dslopes = np.concatenate((axd,azd))

        #Free up some memory
        del axd
        del azd
        del shade
        gc.collect()

    #Create bounds list
    if bounds is None:
        bounds = []
        for i in range(np.shape(ifslopes)[1]):
            bounds.append((smin,smax))

    #Print initial merit function
##    print ampMeritFunction(np.zeros(np.shape(ifslopes)[1]),dslopes,ifslopes)

    #Call optimizer algorithm
    optv = fmin_slsqp(ampMeritFunction,np.zeros(np.shape(ifslopes)[1]),\
                      bounds=bounds,args=(dslopes,ifslopes),\
                      iprint=2,fprime=ampMeritDerivative,iter=200,\
                      acc=1.e-6)

    #Free up optimization memory
    del ifslopes
    del dslopes
    gc.collect()

    #Construct solution image
    ifuncs = pyfits.getdata(ifuncf)
    ishape = np.shape(ifuncs)
    if ifuncs.ndim == 4:
        ifuncs = ifuncs.reshape((ishape[0]*ishape[1],ishape[2],ishape[3]))
    ifuncs = ifuncs.transpose(1,2,0)#*1.e3 #Get to units of microns
    sol = np.dot(ifuncs,optv)

##    pyfits.writeto('Sol.fits',sol)
##    pyfits.writeto('SolVolt.fits',optv)

    return sol,optv

def slopeOptimizer3(ifuncfx,ifuncfy,distortionf,shadef,\
                        azweight=.015,\
                        smin=0.,smax=5.,bounds=None):
    """Assumes IFs and distortion are in slope form in arcsec.
    """
    if dslopes is None:
        #Load in distortion, ifuncs, and shade from filenames
        axif = pyfits.getdata(ifuncfx)
        ayif = pyfits.getdata(ifuncfy)
        shade = pyfits.getdata(shadef)
        
        #Reshape ifuncs into 3D matrix
        ishape = np.shape(axif)
        axif = ifuncx.transpose(1,2,0)
        ayif = ifuncy.transpose(1,2,0)

        #Reshape IFs into 2D matrices
        axif = axif.reshape(((ishape[1]-1)*(ishape[2]-1),ishape[0]))
        azif = azif.reshape(((ishape[1]-1)*(ishape[2]-1),ishape[0]))

        #filter ifuncs array with shade mask
        shade = shade.astype('bool')
        shade = shade.flatten()
        axif = axif[shade,:]
        azif = azif[shade,:]*azweight

        #Merge slopes into 2D optimization matrix
        ifslopes = np.concatenate((axif,azif),axis=0)

        #Create distortion slope arrays
        axd = pyfits.getdata(distortionfx)
        ayd = pyfits.getdata(distortionfy)
        
        #Flatten distortion arrays
        axd = axd.flatten()
        ayd = ayd.flatten()
        axd = axd[shade]
        ayd = ayd[shade]*azweight

        #Create 1D distortion slope target
        dslopes = np.concatenate((axd,ayd))

    #Create bounds list
    if bounds is None:
        bounds = []
        for i in range(np.shape(ifslopes)[1]):
            bounds.append((smin,smax))

    #Print initial merit function
##    print ampMeritFunction(np.zeros(np.shape(ifslopes)[1]),dslopes,ifslopes)

    #Call optimizer algorithm
    optv = fmin_slsqp(ampMeritFunction,np.zeros(np.shape(ifslopes)[1]),\
                      bounds=bounds,args=(dslopes,ifslopes),\
                      iprint=2,fprime=ampMeritDerivative,iter=200,\
                      acc=1.e-6)

    #Construct solution image
    axif = pyfits.getdata(ifuncfx)
    ayif = pyfits.getdata(ifuncfy)
    axif = axif.transpose(1,2,0)
    ayif = ayif.transpose(1,2,0)
    solx = np.dot(axif,optv)
    soly = np.dot(ayif,optv)

##    pyfits.writeto('Sol.fits',sol)
##    pyfits.writeto('SolVolt.fits',optv)

    return solx,soly,optv

#Need a significantly faster merit function
#Option 1: F2PY
#Option 2: Numba jit
#Option 3: CUDA

###Option 2: Numbapro vectorize
##@guvectorize(['void(float32[:],float32[:],float32[:,:],float32[:])'],\
##             '(y),(x),(x,y)->()',\
##             target='parallel')
##def meanSquares(voltages,distortion,ifuncs):
##    #Compute dot product
##    tmp2 = 0.
##    for i in range(np.size(distortion)):
##        tmp = 0.
##        for j in range(np.size(voltages)):
##            tmp += voltages[j]*ifuncs[i,j]
##        tmp2 += (tmp - distortion[i])**2
##    out[:] = tmp2/np.size(distortion)

def slopeMeritFunction(voltages,distortion,ifuncs,dx,azweight):
    """Simple merit function calculator.
    voltages is 1D array of actuator voltages
    distortion is 2D image array
    ifuncs is 3D array, 1st dimension is actuator
    mask has been handled with nans
    New plan is to differentiate and flatten distortion map
    and influence functions prior to optimization
    Then run just like amplitude optimization
    """
    #Numpy way
    resid = distortion - np.dot(ifuncs,voltages)
    ax = np.diff(resid,axis=0)/(dx) #axial slopes
    az = np.diff(resid,axis=1)/(dx) #azimuthal slopes
    
    return np.nansum(ax**2)+np.nansum((az*azweight)**2)

def flatOptimizer(distortion,ifuncs,shade,smin=-5.,smax=5.):
    """Format the arrays to matrices and vectors for
    the merit function. Then call fmin_slsqp to determine
    the optimal voltages to correct the mirror.
    """
    #Reshape ifuncs into 2D matrix
    ishape = np.shape(ifuncs)
    ifuncs = ifuncs.reshape((ishape[0]*ishape[1],ishape[2]*ishape[3]))
    ifuncs = np.transpose(ifuncs)*1000.

    #Flatten shade mask and filter ifuncs array
    shade = shade.flatten()
    shade = shade.astype('bool')
    ifuncs = ifuncs[shade,:]

    #Flatten distortion array and filter with shade mask
    distortion = distortion.flatten()
    distortion = distortion[shade]

    #Create bounds list
    bounds = []
    for i in range(np.shape(ifuncs)[1]):
        bounds.append((smin,smax))

    #Call optimizer algorithm
    optv = fmin_slsqp(ampMeritFunction,np.zeros(np.shape(ifuncs)[1]),bounds=bounds,\
               args=(distortion,ifuncs),iprint=2)

    #Construct solution image
    sol = np.dot(ifuncs,optv)
    img = np.zeros(np.size(shade))
    img[np.where(shade==True)] = sol

    return img,optv

def flatSlopeOptimizer(dslopes=None,ifslopes=None,ifuncf=None,\
                       distortionf=None,shadef=None,\
                       dx=100./150,azweight=.015,\
                       smin=0.,smax=5.,write=False):
    """Format the arrays to matrices and vectors for
    the merit function. Then call fmin_slsqp to determine
    the optimal voltages to correct the mirror.
    All conversion to slope and any frequency filtering must
    occur prior to optimization call.
    """
    #Load in distortion, ifuncs, and shade from file
    #Reshape correctly and then free up memory
    if dslopes is None:
        #Load in distortion, ifuncs, and shade from filenames
        ifuncs = pyfits.getdata(ifuncf)/(dx*1000)*180./np.pi*60**2
        distortion = pyfits.getdata(distortionf)/(dx*1000)*180./np.pi*60**2
        shade = pyfits.getdata(shadef)
        
        #Reshape ifuncs into 3D matrix
        ishape = np.shape(ifuncs)
        ifuncs = ifuncs.reshape((ishape[0]*ishape[1],ishape[2],ishape[3]))
        ifuncs = ifuncs.transpose(1,2,0)*1000. #Get to units of microns
        axif = np.diff(ifuncs,axis=0) #Axial slopes
        axif = axif[:,:-1,:] #Get rid of last column
        azif = np.diff(ifuncs,axis=1) #Azimuthal slopes
        azif = azif[:-1,:,:] #Get rid of last row
        del ifuncs
        gc.collect()

        #Reshape IFs into 2D matrices
        axif = axif.reshape(((ishape[2]-1)*(ishape[3]-1),ishape[0]*ishape[1]))
        azif = azif.reshape(((ishape[2]-1)*(ishape[3]-1),ishape[0]*ishape[1]))

        #filter ifuncs array with shade mask
        shade = shade[:-1,:-1] #Get rid of last row and column
        shade = shade.astype('bool')
        shade = shade.flatten()
        axif = axif[shade,:]
        azif = azif[shade,:]*azweight

        #Merge slopes into 2D optimization matrix
        ifslopes = np.concatenate((axif,azif),axis=0)

        #Free up some memory
        del axif
        del azif
        gc.collect()

        #Create distortion slope arrays
        axd = np.diff(distortion,axis=0)
        axd = axd[:,:-1]
        azd = np.diff(distortion,axis=1)
        azd = azd[:-1,:]
        del distortion
        gc.collect()

        #FLatten distortion arrays
        axd = axd.flatten()
        azd = azd.flatten()
        axd = axd[shade]
        azd = azd[shade]*azweight

        #Create 1D distortion slope target
        dslopes = np.concatenate((axd,azd))

        #Free up some memory
        del axd
        del azd
        del shade
        gc.collect()

    #Create bounds list
    bounds = []
    for i in range(np.shape(ifslopes)[1]):
        bounds.append((smin,smax))

    #Call optimizer algorithm
    optv = fmin_slsqp(ampMeritFunction,np.zeros(np.shape(ifslopes)[1]),\
                      bounds=bounds,args=(dslopes,ifslopes),\
                      iprint=2,fprime=ampMeritDerivative,iter=200,\
                      acc=1.e-6)

    #Free up optimization memory
    del ifslopes
    del dslopes
    gc.collect()

    #Construct solution image
    ifuncs = pyfits.getdata(ifuncf)
    ishape = np.shape(ifuncs)
    ifuncs = ifuncs.reshape((ishape[0]*ishape[1],ishape[2],ishape[3]))
    ifuncs = ifuncs.transpose(1,2,0)*1000. #Get to units of microns
    sol = np.dot(ifuncs,optv)

    if write is True:
        pyfits.writeto('Sol.fits',sol)
        pyfits.writeto('SolVolt.fits',optv)

    return sol,optv

