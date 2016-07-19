import numpy as np
import matplotlib.pyplot as plt

import traces.surfaces as surf
import traces.transformations as tran
import traces.analyses as anal
import traces.sources as sources

import pdb

nSiO2 = 1.4570

def createWavefront(rad,num,coeff,rorder=None,aorder=None):
    """Bounce rays off of Zernike surface. Use flat to
    bring rays to a common plane, leaving the OPD as twice
    the figure error of the Zernike surface.
    """
    #Create set of rays
    rays = sources.circularbeam(rad,num)
    #Reflect to Zernike surface
    surf.zernsurf(rays,coeff,rad,nr=1.,rorder=rorder,aorder=aorder)
    tran.reflect(rays)
    tran.transform(rays,0,0,0,np.pi,0,0)
    surf.flat(rays,nr=1.)
    #Wavefront now has the proper Zernike form, rays pointing in
    #+z direction
    return rays

def traceWedge(rays,t=25.,wang=1.*np.pi/180,pang=45.*np.pi/180):
    """
    Make two copies of rays and trace through a wedged plate.
    Ignore multiple reflections.
    Interpolate one OPD onto the other, take difference
    modulo wavelength
    t = plate thickness (at narrow end)
    ang = wedge angle
    """
    #Make copy
    rays2 = np.copy(rays)

    #Trace first set
    ref1 = [tran.tr.identity_matrix()]*4
    tran.transform(rays,0,0,300.,pang,0,0,coords=ref1)
    surf.flat(rays,nr=1.)
    tran.reflect(rays)
    tran.transform(rays,0,0,0,np.pi/2-pang,0,0,coords=ref1)
    tran.transform(rays,0,0,-300.,0,0,0,coords=ref1)
##    tran.steerY(rays,coords=ref1)
    surf.flat(rays,nr=1.)

    #Trace second set
    ref2 = [tran.tr.identity_matrix()]*4
    pdb.set_trace()
    tran.transform(rays2,0,0,300.,pang,0,0,coords=ref2)
    surf.flat(rays2,nr=1.)
    #Refract into glass and reflect
    tran.refract(rays2,1.,nSiO2)
    tran.transform(rays2,0,0,t,0,wang,0,coords=ref2)
    surf.flat(rays2,nr=nSiO2)
    tran.reflect(rays2)
    #Refract out of glass
##    tran.itransform(rays2,0,0,t,wang,0,0,coords=ref2)
    tran.transform(rays2,0,0,0,0,-wang,0,coords=ref2)
    tran.transform(rays2,0,0,-t,0,0,0,coords=ref2)
    surf.flat(rays2,nr=nSiO2)
    tran.refract(rays2,nSiO2,1.)
    #Go to focal plane
    rays2 = tran.applyT(rays2,ref2,inverse=True)
    rays2 = tran.applyT(rays2,ref1)
    surf.flat(rays2,nr=1.)

    #Both sets of rays at same plane, should have shear and tilt
    #Interpolate OPDs onto common grid
    opd1,dx,dy = anal.interpolateVec(rays,0,200,200,\
                               xr=[rays[1].min(),rays[1].max()],\
                               yr=[rays2[2].min(),rays[2].max()])
    opd2 = anal.interpolateVec(rays2,0,200,200,\
                               xr=[rays[1].min(),rays[1].max()],\
                               yr=[rays2[2].min(),rays[2].max()])[0]

    #Convert to complex phase
    opd1 = opd1/.000635*2*np.pi % (2*np.pi)
    opd2 = opd2/.000635*2*np.pi % (2*np.pi)
    opd1 = np.exp(1j*opd1)
    opd2 = np.exp(1j*opd2)

    #Compute intensity/interferogram
    return np.abs(opd1+opd2)**2
