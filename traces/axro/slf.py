import numpy as np
from numpy import pi,cos,sin,tan,exp,sqrt
import matplotlib.pyplot as plt
import traces.sources as sources
import traces.surfaces as surf
import traces.transformations as tran
import traces.conicsolve as conic
import traces.analyses as anal
import pdb

def sourceToChamber(N,misalign=np.zeros(6)):
    """
    Trace randomly sampled rays from the TruFocus X-ray source
    to the 1.22 m diameter entrance to the test chamber.
    A-B from Jeff K.'s memo is 89.61
    Use oversized sub-apertured annulus, applying translations
    """
    #Define some Wolter parameters
    r1 = conic.primrad(8600.,220.,8400.)
    dphi = 100./220./2
    #Set up subannulus
    rays = sources.subannulus(220.,r1,dphi*1.25,N)
    #Set direction cosines
    srcdist = 89.61e3+(1.5e3-misalign[2])
    raydist = sqrt(srcdist**2+\
                   (rays[1]-misalign[0])**2+\
                   (rays[2]-misalign[1])**2)
    l = (rays[1]-misalign[0])/raydist
    m = (rays[2]-misalign[1])/raydist
    n = -sqrt(1. - l**2 - m**2)
    rays = [rays[0],rays[1],rays[2],rays[3],l,m,n,rays[7],rays[8],rays[9]]
    pdb.set_trace()
    #Go to mirror node and apply rotational misalignment
    tran.transform(rays,220.,0,0,misalign[3],misalign[4],misalign[5])
    tran.transform(rays,-220.,0,0,0,0,0)
    #Place Wolter surfaces
    tran.transform(rays,0,0,-8400.,0,0,0)
    surf.wolterprimary(rays,220.,8400.)
    tran.reflect(rays)
    #Vignette rays not landing in active mirror area
    indz = np.logical_and(rays[3]>8426.,rays[3]<8526.)
    ind = np.logical_and(np.abs(rays[2])<50.,indz)
    rays = tran.vignette(rays,ind=ind)
    #Place secondary
    pdb.set_trace()
    surf.woltersecondary(rays,220.,8400.)
    tran.reflect(rays)
    #Vignette rays not landing in active mirror area
    indz = np.logical_and(rays[3]>8276.,rays[3]<8376.)
    ind = np.logical_and(np.abs(rays[2])<50.,indz)
    rays = tran.vignette(rays,ind=ind)
    #Go back up to intersection plane
    tran.transform(rays,0,0,8400,0,0,0)
    #Reverse misalignments
    tran.itransform(rays,-220.,0,0,0,0,0)
    tran.itransform(rays,0,0,0,misalign[3],misalign[4],misalign[5])
    tran.itransform(rays,220,0,0,0,0,0)
    #Now back in nominal intersection coordinate system
    #Go to focus
    print surf.focusI(rays)
    
    
    return rays

def placeWolterPair(rays,misalign=np.zeros(6)):
    """
    Place the X-ray test mirror pair in the beam path.
    Assume rays are at XY plane with -z mean direction
    Nominal position of intersection plane is 1.5 m
    past chamber entrance with mirror optical axis
    coincident with chamber optical axis.
    Can supply misalignment about X=0,Y=0 in intersection plane.
    """
    #Go to nominal intersection plane
    tran.transform(rays,0,0,-1500.,0,0,0)
    #Apply misalignments
    tran.transform(rays,*misalign)
    #Go to focus and place primary
    tran.transform(rays,0,0,-8400,0,0,0)
    pdb.set_trace()
    surf.wolterprimary(rays,220.,8400.)
    tran.reflect(rays)
    pdb.set_trace()
    #Vignette rays not landing in active mirror area
    indz = np.logical_and(rays[3]>8426.,rays[3]<8526.)
    ind = np.logical_and(np.abs(rays[2])<50.,indz)
    rays = tran.vignette(rays,ind=ind)
    #Place secondary
    pdb.set_trace()
    surf.woltersecondary(rays,220.,8400.)
    tran.reflect(rays)
    #Vignette rays not landing in active mirror area
    indz = np.logical_and(rays[3]>8276.,rays[3]<8376.)
    ind = np.logical_and(np.abs(rays[2])<50.,indz)
    rays = tran.vignette(rays,ind=ind)
    #Go back up to nominal intersection plane
    tran.transform(rays,0,0,8400,0,0,0)
    tran.itransform(rays,*misalign)
    return rays

def findFocus(rays):
    surf.focusX(rays)
    return None
