import numpy as np
import matplotlib.pyplot as plt
import pdb

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
                        detalign=np.zeros(6)):
    """
    Trace rays through the primary mirror and then down to a focus.
    Need to simulate an initial misalignment and then applying
    an optimization algorithm to align primary to beam.
    Merit function should include the random error in spot centroiding
    """
    #Move to primary reference frame - rays 200 mm above node
    tran.transform(rays,0,0,-200.,0,0,0)
    glo = [tran.tr.identity_matrix()]*4
    #Move to mirror tangent point and apply misalignment
    tran.transform(rays,conic.primrad(8450.,220.,8400.),0,50,0,0,0,coords=glo)
    tran.transform(rays,*primalign,coords=glo)
    tran.itransform(rays,conic.primrad(8450.,220.,8400.),0,50,0,0,0,coords=glo)
    tran.transform(rays,0,0,-8400.,0,0,0,coords=glo)
    surf.wolterprimary(rays,220.,8400.)
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
    
    return cen
    
    
