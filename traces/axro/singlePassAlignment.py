import numpy as np
import matplotlib.pyplot as plt
import pdb
import scipy.optimize as opt

import traces.surfaces as surf
import traces.transformations as tran
import traces.analyses as anal
import traces.sources as sources
import traces.conicsolve as conic


def createWavefront(rad,num,coeff,rorder=None,aorder=None,\
                    slitwidth=3.,masknum=15):
    """Bounce rays off of Zernike surface. Use flat to
    bring rays to a common plane, leaving the OPD as twice
    the figure error of the Zernike surface.
    Use subannulus so as not to waste rays in beginning
    of simulation
    Group rays for a given mask slit together using a
    Hartmann vector. Vignette everything else.
    Assume masknum slits 3 mm wide distributed evenly
    over the mirror aperture
    """
    #Create set of rays
    r1 = conic.primrad(8500.,220.,8400.)
    #Loop through Hartmann mask
    maskcenters = np.linspace(-48.5/220.,48.5/220.,masknum)
    for i in range(masknum):
        trays = sources.subannulus(220.,r1,slitwidth/220.,round(num/masknum))
        tran.transform(trays,0,0,0,0,0,maskcenters[i])
        try:
            rays = [np.concatenate([rays[ti],trays[ti]]) for ti in range(10)]
            mask = np.concatenate([mask,np.repeat(i,round(num/masknum))])
        except:
            rays = trays
            mask = np.repeat(i,round(num/masknum))

    tran.transform(rays,220.3,0,0,0,0,0)
    #Reflect to Zernike surface
    surf.zernsurf(rays,coeff,rad,nr=1.,rorder=rorder,aorder=aorder)
    tran.reflect(rays)
    surf.flat(rays,nr=1.)
    tran.transform(rays,-220.3,0,0,0,0,0)
    #Wavefront now has the proper Zernike form, rays pointing in
    #-z direction
    return rays,mask

def traceThroughPrimary(rays,mask,primalign=np.zeros(6),\
                        detalign=np.zeros(6),primCoeffs=None,cenSig=0.):
    """
    Trace rays through the primary mirror and then down to a focus.
    Need to simulate an initial misalignment and then applying
    an optimization algorithm to align primary to beam.
    Merit function should include the random error in spot centroiding
    primCoeffs is a list of coefficients, axial orders, and azimuthal orders
    Use global coordinate systems to determine sign conventions
    """
    #Move to primary reference frame - rays 200 mm above node
    tran.transform(rays,0,0,-200.,0,0,0)
    glo = [tran.tr.identity_matrix()]*4
    #Move to mirror tangent point and apply misalignment
    tran.transform(rays,conic.primrad(8450.,220.,8400.),0,50,0,0,0,coords=glo)
    tran.transform(rays,*primalign,coords=glo)
    tran.itransform(rays,conic.primrad(8450.,220.,8400.),0,50,0,0,0,coords=glo)
    tran.transform(rays,0,0,-8400.,0,0,0,coords=glo)
    #Trace to Wolter surface
    if primCoeffs is None:
        surf.wolterprimary(rays,220.,8400.)
    else:
        surf.primaryLL(rays,220.,8400.,8500.,8400.,100./220.,\
                       *primCoeffs)
    rays = tran.applyT(rays,glo,inverse=True)
    #Rays are now at primary in global coordinate system
    #(origin on optical axis and at nominal node height)
    #Now reflect and trace down to the detector
    tran.reflect(rays)
    tran.transform(rays,0,0,-conic.primfocus(220.,8400.),0,0,0)
    #Apply detector misalignment
    tran.transform(rays,*detalign)
    surf.flat(rays)
    #Pick out spot centroids
    cen = [anal.centroid(rays,weights=mask==i) for i in range(mask[-1]+1)]
    cen = np.transpose(np.array(cen))
    #Add centroiding error
    if cenSig > 0:
        cen = cen + np.random.normal(scale=cenSig,size=np.shape(cen))

    return cen
    
def primaryTrace(rad,num,coeff,primalign=np.zeros(6),detalign=np.zeros(6)):
    """
    Function to create source rays and trace through the primary
    to the focus detector.
    """
    #Trace
    rays,mask = createWavefront(rad,num,coeff)
    cen = traceThroughPrimary(rays,mask,primalign=primalign,detalign=detalign)
    #Return deviation from centroid of centroid
    cx = np.mean(cen[0])
    cy = np.mean(cen[1])
    return np.sqrt(np.mean((cx-cen[0])**2+(cy-cen[1])**2))

def alignPrimary(rad,num,coeff,initial=np.zeros(6),detalign=np.zeros(6)):
    """
    Simulate the process of aligning the primary mirror to the
    plane wave. Assume a starting misalignment of the primary as
    input to the optimizer, and then allow the angular degrees
    of freedom to vary.
    Most of this can be done using scipy.optimize.minimize
    Question of which algorithm to use, likely Nelder-Mead
    """
    fun = lambda p: primaryTrace(rad,num,coeff,\
                                 primalign=[initial[0],initial[1],initial[2],\
                                            p[0],p[1],initial[-1]],\
                                 detalign=detalign)
    
    res = opt.minimize(fun,initial[3:5],method='Nelder-Mead')

    return res

def traceThroughPair(rays,mask,primalign=np.zeros(6),\
                        detalign=np.zeros(6),primCoeffs=None,cenSig=0.):
    """
    Trace rays through the primary mirror and then down to a focus.
    Need to simulate an initial misalignment and then applying
    an optimization algorithm to align primary to beam.
    Merit function should include the random error in spot centroiding
    primCoeffs is a list of coefficients, axial orders, and azimuthal orders
    Use global coordinate systems to determine sign conventions
    """
    #Move to primary reference frame - rays 200 mm above node
    tran.transform(rays,0,0,-200.,0,0,0)
    glo = [tran.tr.identity_matrix()]*4
    #Move to mirror tangent point and apply misalignment
    tran.transform(rays,conic.primrad(8450.,220.,8400.),0,50,0,0,0,coords=glo)
    tran.transform(rays,*primalign,coords=glo)
    tran.itransform(rays,conic.primrad(8450.,220.,8400.),0,50,0,0,0,coords=glo)
    tran.transform(rays,0,0,-8400.,0,0,0,coords=glo)
    #Trace to Wolter surface
    if primCoeffs is None:
        surf.wolterprimary(rays,220.,8400.)
    else:
        surf.primaryLL(rays,220.,8400.,8500.,8400.,100./220.,\
                       *primCoeffs)
    rays = tran.applyT(rays,glo,inverse=True)
    #Rays are now at primary in global coordinate system
    #(origin on optical axis and at nominal node height)
    #Now reflect and move on to secondary
    tran.reflect(rays)
    glo = [tran.tr.identity_matrix()]*4
    tran.transform(rays,conic.secrad(8350.,220.,8400.),0,-50.,0,0,0,coords=glo)
    tran.transform(rays,*secalign,coords=glo)
    tran.itransform(rays,conic.secrad(8350.,220.,8400.),0,-50.,0,0,0,coords=glo)
    tran.transform(rays,0,0,-8400.,0,0,0,coords=glo)
    #Trace to secondary surface
    if secCoeffs is None:
        surf.woltersecondary(rays,220.,8400.)
    else:
        surf.secondaryLL(rays,220.,8400.,8400.,8300.,100./220.,\
                         *secCoeffs)
    rays = tran.applyT(rays,glo,inverse=True)
    #Rays are now at secondary in global coordinate system
    #(origin on optical axis and at nominal node height)
    #Now reflect and go to detector
    tran.reflect(rays)
