import numpy as np
import matplotlib.pyplot as plt
from numpy import exp,cos,sin,pi,sqrt
import traces.PyTrace as PT

nbk7 = 1.5150885
nsf2 = 1.6437928

def collimation(rays,reverse=False):
    """Assuming reference frame is +Z toward optical axis and
    XY plane tangent with first surface, trace through the
    collimation lens
    R1=9324.
    R2=1124.
    Thickness=20.
    """
    R1 = 9324.
    R2 = 1124.
    t = 20.
    PT.transform(rays,0,0,R1,0,0,0) #Center of first surface
    PT.sphere(rays,9324.,nr=1.) #Trace to surface
    PT.refract(rays,1.,nbk7) #Refract into material
    PT.transform(rays,0,0,-R1+t-R2,0,0,0) #Center of second surface
    PT.sphere(rays,R2,nr=nbk7) #Trace to second surface
    PT.refract(rays,nbk7,1.) #Refract back to air
    PT.transform(rays,0,0,R2,0,0,0) #Set reference frame tangent to last surface
    PT.flat(rays,nr=1.) #Trace to tangent plane
    return
