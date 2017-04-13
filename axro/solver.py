import numpy as np
#from numbapro.cudalib import cublas
import pdb,os
from scipy.optimize import fmin_slsqp,least_squares
#from numbapro import guvectorize
#from axro.merit import merit
import astropy.io.fits as pyfits
##import gc
import scipy.interpolate as interp
import utilities.transformations as tr
import traces.conicsolve as conic
import time

#import cvxpy

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

def prepareIFs(ifs,dx=None,azweight=.015):
    """
    Put IF arrays in format required by optimizer.
    If dx is not None, apply derivative.
    """
    #Apply derivative if necessary
    #First element of result is axial derivative
    if dx is not None:
        ifs = np.array(np.gradient(ifs,*dx,axis=(1,2)))*180/np.pi*60.**2 / 1000.
        ifs[1] = ifs[1]*azweight
        ifs = ifs.transpose(1,0,2,3)
        sha = np.shape(ifs)
        for i in range(sha[0]):
            for j in range(sha[1]):
                ifs[i,j] = ifs[i,j] - np.nanmean(ifs[i,j])
        ifs = ifs.reshape((sha[0],sha[1]*sha[2]*sha[3]))
    else:
        #ifs = ifs.transpose(1,2,0)
        sha = np.shape(ifs)
        for i in range(sha[0]):
            ifs[i] = ifs[i] - np.nanmean(ifs[i])
        ifs = ifs.reshape((sha[0],sha[1]*sha[2]))

    return np.transpose(ifs)

def prepareDist(d,dx=None,azweight=.015):
    """
    Put distortion array in format required by optimizer.
    If dx is not None, apply derivative.
    Can also be run on shademasks
    """
    #Apply derivative if necessary
    #First element of result is axial derivative
    if dx is not None:
        d = np.array(np.gradient(d,*dx))*180/np.pi*60.**2 / 1000.
        d[0] = d[0] - np.nanmean(d[0])
        d[1] = d[1] - np.nanmean(d[1])
        d[1] = d[1]*azweight

    return d.flatten()

def optimizer(distortion,ifs,shade,smin=0.,smax=5.,bounds=None,compare=False):
    """
    Cleaner implementation of optimizer. ifs and distortion should
    already be in whatever form (amplitude or slope) desired.
    IFs should have had prepareIFs already run on them.
    Units should be identical between the two.
    """
    #Load in data
    if type(distortion)==str:
        distortion = pyfits.getdata(distortion)
    if type(ifs)==str:
        ifs = pyfits.getdata(ifs)
    if type(shade)==str:
        shade = pyfits.getdata(shade)

    #Remove shademask
    ifs = ifs[shade==1]
    distortion = distortion[shade==1]

    #Remove nans
    ind = ~np.isnan(distortion)
    ifs = ifs[ind]
    distortion = distortion[ind]

    if compare is True:
        #Output arrays as fits to be compared using MATLAB
        return ifs,distortion
    

    #Handle bounds
    if bounds is None:
        bounds = []
        for i in range(np.shape(ifs)[1]):
            bounds.append((smin,smax))

    #Call optimizer algorithm
    tstart = time.time()
    optv = fmin_slsqp(ampMeritFunction,np.zeros(np.shape(ifs)[1]),\
                      bounds=bounds,args=(distortion,ifs),\
                      iprint=1,fprime=ampMeritDerivative,iter=200,\
                      acc=1.e-6)
    print time.time()-tstart

    return optv

def correctDistortion(dist,ifs,shade,dx=None,azweight=.015,smax=5.,\
                      bounds=None,compare=False):
    """
    Wrapper function to apply and evaluate a correction
    on distortion data.
    Distortion and IFs are assumed to already be on the
    same grid size.
    dx should be in mm, dist and ifs should be in microns
    """
    #Make sure shapes are correct
    if not (np.shape(dist)==np.shape(ifs[0])==np.shape(shade)):
        print 'Unequal shapes!'
        return None

    #Prepare arrays
    distp = prepareDist(dist,dx=dx,azweight=azweight)
    ifsp = prepareIFs(ifs,dx=dx,azweight=azweight)
    shadep = prepareDist(shade)

    #Run optimizer
    res = optimizer(-distp,ifsp,shadep,smax=smax,bounds=bounds,compare=compare)

##    #Try with CVXPY
##    tstart = time.time()
##    m,n = np.shape(ifsp)
##    v = cvxpy.Variable(n)
##    objective = cvxpy.Minimize(cvxpy.sum_entries(cvxpy.square(ifsp*v + distp)))
##    constraints = [0 <= v, v<=5.]
##    prob = cvxpy.Problem(objective,constraints)
##
##    prob.solve()
##
##    print prob.value
##    print v.value
##    print 'CVXPY: ' + str(time.time()-tstart)

##    #Reconstruct solution and initial slope error
##    gdist = np.gradient(dist,*dx)
##    ifs = ifs.transpose(1,2,0)
##    correction = np.dot(ifs,res)
##    gcor = np.gradient(correction,*dx)
##    resid = gdist[0] + gcor[0]

    return res

def testcase():
    ifs = pyfits.getdata('/home/rallured/Dropbox/AXRO/Simulation/'
                         'MatlabOptimization/170320_TestIFs.fits')
    distort = pyfits.getdata('/home/rallured/Dropbox/AXRO/Simulation'
                             '/MatlabOptimization/170320_TestDis.fits')
    res = optimizer(distort,ifs,np.ones(len(distort)))
    return res
    

solvedir = '/home/rallured/solve_pzt/'
def setupAlexey(iffile,dx,shadefile=solvedir+'shademasks/Pyshade.fits',\
                disfile=solvedir+'distortions/Pydist.fits',\
                pzt_max_strain=5.,pzt_min_strain=0.,minmode='slope',\
                azim_slope_weight=.015):
    """
    Set up a .par file for solve_pzt
    """
    #Write python.par
    f = open('/home/rallured/solve_pzt/parfiles/python.par','w')
    f.write('influence_function='+iffile+'\n')
    f.write('influence_function_symstorage=no\n'
            'influence_function_reflect_y=no\n'
            'influence_function_reflect_x=no\n'
            'influence_function_units=mum\n'
            'pixel_size='+str(dx)+' mm\n'
            'pzt_max_strain='+str(pzt_max_strain)+'\n'
            'pzt_min_strain='+str(pzt_min_strain)+'\n'
            'neg_in_shade=no\n'
            'shade_map='+str(shadefile)+'\n'
            'fit_tilt=yes\n'
            'fit_offset=no\n'
            'minmode='+str(minmode)+'\n'
            'azim_slope_weight='+str(azim_slope_weight)+'\n'
            'distort_map='+str(disfile)+'\n')
    f.close()

    return None

def preMath():
    """
    Execute python.par and save pre_math
    """
    #Go to solve_pzt directory and run solver
    pwd = os.getcwd()
    os.chdir('/home/rallured/solve_pzt/parfiles')
    os.system('~/solve_pzt/smartx/bin/solve_pzt @python.par '
              'save_premath=Pymath.dat -o Py')

    #Go back to working directory
    os.chdir(pwd)

    return

def pyExecute(distortion,disfile=solvedir+'distortions/Pydist.fits'):
    """
    Write distortion image into disfile and execute
    solve_pzt, fetching results
    """
    #Write distortion map
    pyfits.writeto(disfile,distortion,clobber=True)
    
    #Run solver
    pwd = os.getcwd()
    os.chdir('/home/rallured/solve_pzt/parfiles')
    os.system('~/solve_pzt/smartx/bin/solve_pzt @python.par '
              'premath=Pymath.dat -o Py')
    #Get results
    cor = pyfits.getdata('Py-sol.fits')
    resid = pyfits.getdata('Py-resid.fits')
    volt = pyfits.getdata('Py-adj-map.fits')
    
    #Go back to working directory             
    os.chdir(pwd)

    return cor,resid,volt

def convertFEAInfluence(filename,Nx,Ny,method='cubic',\
                        cylcoords=True):
    """Read in Vanessa's CSV file for AXRO mirror
    Mirror no longer assumed to be cylinder.
    Need to regrid initial and perturbed nodes onto regular grid,
    then compute radial difference.
    """
    #Load FEA data
    d = np.transpose(np.genfromtxt(filename,skip_header=1,delimiter=','))

    if cylcoords is True:
        r0 = d[2]*1e3
        rm = np.mean(r0)
        t0 = d[3]*np.pi/180. * rm #Convert to arc length in mm
        z0 = d[4]*1e3
        #r0 = np.repeat(220.497,len(t0))
        r = r0 + d[5]*1e3
        t = (d[3] + d[6])*np.pi/180. * rm #Convert to arc length in mm
        z = z0 + d[7]*1e3
    else:
        x0 = d[2]*1e3
        y0 = d[3]*1e3
        z0 = d[4]*1e3
        x = x0 + d[5]*1e3
        y = y0 + d[6]*1e3
        z = z0 + d[7]*1e3
        #Convert to cylindrical
        t0 = np.arctan2(x0,-z0)*220.497 #Convert to arc length in mm
        r0 = np.sqrt(x0**2+z0**2)
        z0 = y0
        t = np.arctan2(x,-z)*220.497
        r = np.sqrt(x**2+z**2)
        z = y

    #Construct regular grid
    gy = np.linspace(z0.min(),z0.max(),Nx+2)
    gx = np.linspace(t0.min(),t0.max(),Ny+2)
    gx,gy = np.meshgrid(gx,gy)

    #Run interpolation
    g0 = interp.griddata((z0,t0),r0,(gy,gx),method=method)
    g0[np.isnan(g0)] = 0.
    g = interp.griddata((z,t),r,(gy,gx),method=method)
    g[np.isnan(g)] = 0.

    print filename + ' done'
    
    return -(g0[1:-1,1:-1]-g[1:-1,1:-1]),g0[1:-1,1:-1],g[1:-1,1:-1]

def createShadePerimeter(sh,axialFraction=0.,azFraction=0.):
    """
    Create a shademask where a fraction of the axial and
    azimuthal perimeter is blocked. Fraction is total fraction
    blocked.
    sh is shape tuple e.g. (200,200)
    """
    arr = np.zeros(sh)
    axIndex = round(sh[0]*axialFraction/2)
    azIndex = round(sh[1]*azFraction/2)
    arr[axIndex:-axIndex,azIndex:-azIndex] = 1.
    return arr
