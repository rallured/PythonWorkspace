import numpy as np
import matplotlib.pyplot as plt

import traces.surfaces as surf
import traces.transformations as tran
import traces.analyses as anal
import traces.sources as sources
import traces.conicsolve as conic

def traceZeta(pmin,R0=220.,Z0=1e4,psi=1.,offaxis=0.):
    """
    Set initial aperture based on height of bottom of mirror
    above node (pmin). Assume 100 mm long mirrors.
    Then, trace to secondary and mark smin and smax for
    where the rays strike.
    Then trace out off axis field positions and determine
    RMS and HPD vs angle.
    """
    #Set up aperture
    R0 = conic.primrad(pmin,R0,Z0)
    R1 = conic.primrad(pmin+100.,R0,Z0)
    rays = sources.annulus(R0,R1,1e4)
    tran.transform(rays,0,0,-Z0,0,0,0)

    #Trace to primary and add off-axis angle
    surf.wsPrimary(rays,R0,Z0,psi)
    rays[5] = np.sin(offaxis)
    rays[6] = -np.sqrt(1.-rays[5]**2)
    tran.reflect(rays)

    #Trace to secondary
    surf.wsSecondary(rays,R0,Z0,psi)
    tran.reflect(rays)
    smax = np.nanmax(rays[3])
    smin = np.nanmin(rays[3])

    #Go to focus
    f = surf.focusI(rays)

    #Compute merit functions
    hpd = anal.hpd(rays)
    rms = anal.rmsCentroid(rays)

    return smin,smax,f,hpd,rms
