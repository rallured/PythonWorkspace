import numpy as np
import matplotlib.pyplot as plt
import traces.PyTrace as PT
import pdb

#Set up incident beam trace and determine sensitivity to beam
#impact location.
#Trace nominal geometry (function of incidence angle) and
#record location of diffracted and reflected spots
#Perturb location and angle of beam impact and record
#change in spot locations

#What don't you know about geometry?
#Location and angle of source -> location and angle of beam impact
#Location of detectors (angle is very small contributor)

#Write raytrace with beam impact orientation and
#grating orientation as free parameters

def alignTrace(inc,impact,grating,detector,order=0):
    """Traces UV laser rays to grating. Beam impact misalignment
    is handled with a single coordinate transformations right after
    source definition. Grating orientation is handled with
    symmetric coordinate transformations.
    inc - nominal beam glancing angle, must be less than
        50.39 deg for 262 nm light
    impact - 6 element array giving beam impact transform
    grating - 6 element array giving grating misalignment
    """
    #Set up source with single ray, diffraction plane
    #is XZ, glancing angle from XY plane, ray starts out
    #pointing +x and -z
    PT.pointsource(0.,1)
    PT.transform(0,0,0,0,-np.pi/2-inc,0)
    #Perform beam impact misalignment transform, rotation first
    PT.transform(*np.concatenate(((0,0,0),impact[3:])))
    PT.transform(*np.concatenate((impact[:3],(0,0,0))))
    #Perform grating misalignment
    PT.transform(*grating)
    #Linear grating
    PT.flat()
    PT.reflect()
    PT.grat(160.,order,262.)
    #Reverse misalignment transformation
    PT.itransform(*grating)
    #Go to detector depending on order
    if order is not 0:
        PT.transform(-800.,0,0,0,0,0)
    else:
        PT.transform(800.,0,0,0,0,0)
    #Trace to detector
    PT.transform(0,0,0,0,-np.pi/2,0)
    PT.transform(*detector)
    PT.flat()
    #Return ray position
    return PT.x,PT.y

def computeYaw(inc,impact,grating,detr,detd):
    """Traces both orders and computes yaw of grating
    Uses alignTrace
    Returns yaw angle assuming 
    """
    xr,yr = alignTrace(inc,impact,grating,detr,order=0)
    betar = 800./np.sqrt(800.**2+xr**2+yr**2)
    alphar = -yr/np.sqrt(800.**2+xr**2+yr**2)
    xd,yd = alignTrace(inc,impact,grating,detd,order=1)
    betad = -800./np.sqrt(800.**2+xd**2+yd**2)
    alphad = -yd/np.sqrt(800.**2+xd**2+yd**2)
    
    return np.arctan((alphad-alphar)/(betar-betad))*180./np.pi*60.**2


def dofSensitivity(inc,alignvector,obj='beam',dof=0):
    """Compute x,y positions of reflected and diffracted spots
    return as a function of alignvector
    """
    #Initialize misalignment vectors
    grating = np.zeros(6)
    impact = np.zeros(6)

    #Initialize output vectors
    xr = np.zeros(np.size(alignvector))
    yr = np.copy(xr)
    xd = np.copy(xr)
    yd = np.copy(xr)

    #Perform raytraces in loop
    for a in alignvector:
        #Adjust misalignments
        if obj is 'beam':
            impact[dof] = a
        else:
            grating[dof] = a

        #Perform trace and set appropriate output elements
        i = a==alignvector
        x,y = alignTrace(inc,impact,grating,order=0)
        xr[i] = x
        yr[i] = y
        x,y = alignTrace(inc,impact,grating,order=1)
        xd[i] = x
        yd[i] = y

    #Return
    return xr,yr,xd,yd

