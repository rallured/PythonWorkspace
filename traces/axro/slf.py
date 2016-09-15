import numpy as np
from numpy import pi,cos,sin,tan,exp,sqrt
import matplotlib.pyplot as plt
import traces.sources as sources
import traces.surfaces as surf
import traces.transformations as tran
import traces.conicsolve as conic
import traces.analyses as anal
import pdb
import scipy.signal as sig

def singleOptic(N,misalign=np.zeros(6)):
    """Trace single primary mirror from SLF finite
    source distance.
    """
    #Define some Wolter parameters
    r1 = conic.primrad(8600.,220.,8400.)
    dphi = 100./220./2
    #Set up subannulus
    rays = sources.subannulus(220.,r1,dphi*1.25,N)
##    #Set direction cosines
##    srcdist = 89.61e3+(1.5e3-misalign[2])
##    raydist = sqrt(srcdist**2+\
##                   (rays[1]-misalign[0])**2+\
##                   (rays[2]-misalign[1])**2)
##    l = (rays[1]-misalign[0])/raydist
##    m = (rays[2]-misalign[1])/raydist
##    n = -sqrt(1. - l**2 - m**2)
##    rays = [rays[0],rays[1],rays[2],rays[3],l,m,n,rays[7],rays[8],rays[9]]
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
    #Go to focus
    surf.flat(rays)
    f = surf.focusI(rays)-8400.

    return rays#anal.hpd(rays)/abs(f)*180/pi*60**2,abs(f)

def singleOptic2(n,misalign=np.zeros(6),srcdist=89.61e3+1.5e3,az=50.,\
                 returnRays=False,f=None):
    """Alternative SLF finite source trace"""
    #Establish subannulus of rays
    rays = sources.subannulus(220.,221.,az*1.2/220.,n,zhat=-1.)
    #Transform to node position
    tran.transform(rays,220,0,0,0,0,0)
    #Set up finite source distance
    raydist = sqrt(srcdist**2+rays[1]**2+rays[2]**2)
    l = rays[1]/raydist
    m = rays[2]/raydist
    n = -sqrt(1.-l**2-m**2)
    rays = [raydist,rays[1],rays[2],rays[3],l,m,n,rays[7],rays[8],rays[9]]
    #Align perfectly to beam
    tran.steerX(rays)
    #Apply misalignment
    tran.transform(rays,*misalign)
    #Place mirror
    surf.wolterprimarynode(rays,220,8400.)
    #Vignette rays not landing in active mirror area
    indz = np.logical_and(rays[3]>26.,rays[3]<126.)
    ind = np.logical_and(np.abs(rays[2])<az/2.,indz)
    rays = tran.vignette(rays,ind=ind)
    #Reverse misalignment
    tran.itransform(rays,*misalign)
    #Reflect and go to surface
    tran.reflect(rays)
    if f is None:
        f = surf.focusI(rays)
    else:
        tran.transform(rays,0,0,f,0,0,0)
        surf.flat(rays)
    #Get centroid
    cx,cy = anal.centroid(rays)

    if returnRays is True:
        return rays
    
    return anal.hpd(rays)/abs(f)*180/pi*60**2,f,cx

def examineAzimuthalStripSize(strip=np.linspace(2.5,15.,5)):
    pitch = np.linspace(0,10*.3e-3,200)

    foc = []
    perf = []
    lat = []
    angle = []
    yawbound = []
    pbound = []
    dbound = []
    for s in strip:
        #Compute performance as function of pitch
        res = np.transpose(\
            np.array(\
                [singleOptic2(10000,misalign=[0,0,0,0,t,0],\
                              az=s) for t in pitch]))
        #Find optimal performance
        per = sig.savgol_filter(res[0],11,3)
        ind = np.argmin(per)
        angle.append(pitch[ind])
        perf.append(res[0][ind])
        foc.append(res[1][ind])
        lat.append(res[2][ind])
        #Find alignment sensitivities
        drange = np.linspace(-500.,500.,101)
        angrange = np.linspace(-.3e-3,.3e-3,101)
        yawres = np.transpose(np.array([singleOptic2(1000,misalign=\
                                                  [0,0,0,y,angle[-1],0],\
                                                  f=foc[-1],az=s) \
                                     for y in angrange]))
        plt.figure('yaw')
        plt.plot(angrange*180/pi*60**2,sig.savgol_filter(yawres[0],11,3),\
                 label=str(s))
        yawbound.append(abs(angrange[np.argmin(abs(yawres[0]-1.))]))
        pres = np.transpose(np.array([singleOptic2(1000,misalign=\
                                                  [0,0,0,0,p+angle[-1],0],\
                                                  f=foc[-1],az=s) \
                                     for p in angrange]))
        plt.figure('pitch')
        plt.plot(angrange*180/pi*60**2,sig.savgol_filter(pres[0],11,3),\
                 label=str(s))
        pbound.append(abs(angrange[np.argmin(abs(pres[0]-1))]))
        dres = np.transpose(np.array([singleOptic2(1000,misalign=\
                                                  [0,0,0,0,angle[-1],0],\
                                                  f=foc[-1]+d,az=s) \
                                     for d in drange]))
        plt.figure('d')
        plt.plot(drange,sig.savgol_filter(dres[0],11,3),label=str(s))
        dbound.append(abs(drange[np.argmin(abs(dres[0]-1))]))

    return np.array([perf,foc,lat,angle,yawbound,pbound,dbound])
        
        

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
    f = -9253.3858232
    tran.transform(rays,0,0,f,0,0,0)
    surf.flat(rays)
    
    return rays#anal.hpd(rays)/abs(f)*60**2*180/pi

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
