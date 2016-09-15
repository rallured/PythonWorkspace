import numpy as np
from numpy import cos,sin,sqrt,pi
import matplotlib.pyplot as plt
import pdb

import traces.sources as source
import traces.surfaces as surf
import traces.analyses as anal
import traces.transformations as tran

def traceSPO(N,rin=700.,rout=737.,azwidth=66.,srcdist=89.61e3+1.5e3):
    """
    Trace a set of rays through an SPO module using a
    finite source distance. Ignore vignetting, we are
    only interested in aberrations.
    Set up a subannulus, apply finite source effect,
    and then simply translate to inner SPO radius and
    trace outward.
    Let x be the radial direction, y the azimuthal
    """
    #Establish subannulus of rays
    rays = source.subannulus(rin,rout,azwidth/rin,N,zhat=-1.)
    #Transform to node position
    mx = np.mean([rin,rout])
    tran.transform(rays,mx,0,0,0,0,0)
    #Set up finite source distance
    raydist = np.sqrt(srcdist**2+rays[1]**2+rays[2]**2)
    l = rays[1]/raydist
    m = rays[2]/raydist
    n = -np.sqrt(1.-l**2-m**2)
    rays = [raydist,rays[1],rays[2],rays[3],l,m,n,rays[7],rays[8],rays[9]]
    #Align perfectly to beam
    tran.steerX(rays)
    
    #Move to SPO optical axis and trace through shells
    tran.transform(rays,-mx,0,0,0,0,0)
    R = np.arange(rin,rout+.605,.605)
    rad = np.sqrt(rays[1]**2+rays[2]**2)
    for r in R:
        #Collect relevant rays
        ind = np.logical_and(rad>r,rad<r+.605)
        if np.sum(ind)==0:
            continue
        #Trace them through system
        surf.spoPrimary(rays,r,12e3,ind=ind)
        tran.reflect(rays,ind=ind)
        surf.spoSecondary(rays,r,12e3,ind=ind)
        tran.reflect(rays,ind=ind)
    #Rays are now at secondary surfaces, 
    return rays

def traceOPG(rays,hubdist=11832.911,yaw=0.,order=1,wave=1.,ang=2.5/11832.911):
    """
    Trace the OPG module. Probably ignore vignetting again.
    Place perfect OPG surfaces at the correct angular distance
    to make this a reasonable approximation.
    Assume reference frame is in center of module with -z
    pointing toward hub - achieved with steerX/steerY and
    rotate inc before this function call
    Create vector to keep track of which grating each ray
    diffracts from. Separate LSFs can be identified using
    this vector.
    """
    #Establish starting coordinate system
    coords = [tran.tr.identity_matrix()]*4
    #Get -x pointing to hub
    #Question whether to rotate about z to swap x and y
    tran.transform(rays,0,0,0,0,0,-pi/2,coords=coords)
    tran.transform(rays,0,0,0,pi/2,0,0,coords=coords)
    #Go to hub, then rotate to extreme grating surface
    tran.transform(rays,0,0,0,0,0,yaw,coords=coords) #possible blaze
    tran.transform(rays,0,-11832.911,0,0,0,0,coords=coords)
    tran.transform(rays,0,0,0,-ang*7,0,0,coords=coords) #minus sign ambiguity
    #Loop through gratings, tracing rays
    left = np.repeat(True,len(rays[1]))
    record = np.zeros(len(rays[1]))
    for i in range(15):
        print i
        #If no rays left, we are done
        if np.sum(left) == 0:
            continue
        #Rays with small incidence angle removed
        indg = np.abs(np.arcsin(rays[6])) > .001
        ind = np.logical_and(left,indg)
        if np.sum(ind)==0:
            tran.transform(rays,0,0,0,ang,0,0,coords=coords)
            continue
        #Trace rays to surface
        surf.flat(rays,ind=ind)
        #Identify relevant rays
        ind = np.logical_and(rays[2]>11832.911-96./2,rays[2]<11832.911+96./2)
        ind = np.logical_and(ind,left)
        #Remove these rays from the set that remain
        left = np.logical_and(left,np.invert(ind))
        if np.sum(ind)==0:
            tran.transform(rays,0,0,0,ang,0,0,coords=coords)
            continue
        #Record which grating these rays diffracted from
        record[ind] = i+1
        #Diffract this set of rays
        tran.reflect(rays,ind=ind)
        tran.transform(rays,0,11832.911-hubdist,0,0,0,0,coords=coords)
        tran.radgrat(rays,160./hubdist,order,wave,ind=ind)
        tran.transform(rays,0,hubdist-11832.911,0,0,0,0,coords=coords)
        #Rotate to next grating
        tran.transform(rays,0,0,0,ang,0,0,coords=coords)
    #Go back to original coordinate system
    rays = tran.applyT(rays,coords,inverse=True)

    return rays,record

def test(N,rin=700.,rout=737.,azwidth=66.,srcdist=89.61e3+1.5e3,\
         hubdist=11832.911,yaw=0.,wave=6.,order=1,\
         opgalign=[0,0,0,0,0,0],f=None,\
         rrays=False,glob=False):
    """
    Trace through the SPO module, then place the OPG module
    at its nominal position, allowing for misalignments about the
    center of the OPG module. The module tolerances can be
    investigated by a coordinate transformation around the
    OPG module placement.
    """
    #Trace through SPO module
    rays = traceSPO(N,rin=rin,rout=rout,azwidth=azwidth,srcdist=srcdist)
    coords = [tran.tr.identity_matrix()]*4

    #Find the nominal OPG module location using formalism
    #from Flanagan's SPIE paper
    #Go to focus, steer out X and Y, then go up a distance
    #defined using Flangan formula, this should leave you
    #at the center of the beam, therefore the center of the
    #OPG module
    tran.steerX(rays,coords=coords)
    tran.steerY(rays,coords=coords)
    f0 = surf.focusI(rays,coords=coords)
    tran.transform(rays,np.mean(rays[1]),np.mean(rays[2]),0,0,0,0,\
                   coords=coords)
    tran.transform(rays,0,0,0,0,pi,0,coords=coords)
    tran.transform(rays,0,0,11832.911*np.cos(1.5*np.pi/180),0,0,0,coords=coords)
    tran.transform(rays,0,0,0,0,1.5*np.pi/180,0,coords=coords)
    surf.flat(rays)
    #Now at center of central grating, with -z pointing toward hub
    tran.transform(rays,*opgalign,coords=coords)
    rays,record = traceOPG(rays,hubdist=hubdist,yaw=yaw,wave=wave,order=order)
    tran.itransform(rays,*opgalign,coords=coords)
    #Should be at same reference frame, with rays now diffracted
    if np.sum(record)==0:
        pdb.set_trace()
    rays = tran.vignette(rays,ind=record>0)
    record = record[record>0]

    #Trace to detector and determine LSF
    rays = tran.applyT(rays,coords,inverse=True)
    #surf.focusI(rays)
    if f is not None:
        try:
            tran.transform(rays,0,0,-f,0,0,0)
            surf.flat(rays)
        except:
            pdb.set_trace()

    if rrays is True:
        if glob is True:
            tran.transform(rays,0,0,f,0,0,0)
        return rays,record

    #Return LSF in arcseconds
    return anal.hpdY(rays)/12e3*180/pi*60**2
