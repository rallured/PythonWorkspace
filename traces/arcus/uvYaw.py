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

def alignTrace(inc,impact,grating,order=0):
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
        PT.transform(-400.,0,0,0,0,0)
    else:
        PT.transform(400.,0,0,0,0,0)
    #Trace to detector
    PT.transform(0,0,0,0,np.pi/2,0)
    PT.flat()
    #Return ray position
    return PT.x,PT.y

def dofSensitivity(inc):
    return

