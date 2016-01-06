import numpy as np
import matplotlib.pyplot as plt
from numpy import exp,cos,sin,pi,sqrt
import traces.transformations as tran
import traces.sources as sources
import traces.lenses as lenses
import traces.analyses as anal
import traces.surfaces as surf
import pdb

nbk7 = 1.5150885
nsf2 = 1.6437928
nsil = 1.45701735

#Get CGH coefficients
cghcoeff = np.genfromtxt('/home/rallured/Dropbox/'
                         'AXRO/Metrology/200mmCGH.txt')[2:]
#Empirically determined line focus distance
line = 199.997297
#Empirically determined focus from field lens
foc = 264.92386077410214#275.01818439272461
#Empriically determined position of cylindrical field lens
cylz = 120.70707070707071

def focusCyl(cylz,div=.1*pi/180):
    raylist = cylindricalSource(N=1000,div=div)
    rmsy = [backToWFS(rays,cylz) for rays in raylist]
    return np.mean(rmsy)


def rayBundle(N,div,az,height):
    """
    Set up a diverging ray bundle on the 220 mm cylinder.
    """
    #Establish rays
    rays = sources.pointsource(div,N)
    #Go to cylindrical axis
    tran.transform(rays,0,0,220.,0,0,0)
    #Apply height offset
    tran.transform(rays,-height,0,0,0,0,0)
    #Apply azimuthal offset
    tran.transform(rays,0,0,0,-az,0,0)
    #Go back to tangent plane
    tran.transform(rays,0,0,-220.,0,0,0)
    surf.flat(rays,nr=1.)
    return rays

def cylindricalSource(height=np.linspace(-50.,50.,3),\
                      az=np.linspace(-pi/12,pi/12,5),\
                      N=100,div=pi/180.):
    """
    Set up source rays for focus testing. Set them up on the
    220 mm nominal cylindrical surface using coordinate transforms.
    Use extreme heights and azimuths as well as central sources.
    Use N rays per position, with div divergence
    """
    raylist = [rayBundle(N,div,a,h) for a in az for h in height]
    return raylist

def backToWFS(rays):
    """
    Trace rays from nominal test optic tangent plane back to WFS plane.
    This function can also be used with a point source to determine the
    Optimal focus positions of the field lenses.
    +z points toward CGH.
    """
    #Back to CGH
    tran.transform(rays,0,0,220+line,0,0,0)
    surf.flat(rays,nr=1.)
    #Trace back through CGH
    tran.transform(rays,0,0,0,0,1.*pi/180,0)
    tran.transform(rays,0,0,0,1.*pi/180,0,0)
    surf.flat(rays,nr=1.)
    surf.zernphase(rays,cghcoeff,80.,632.82e-6)
    tran.refract(rays,1.,nsil)
    tran.transform(rays,0,0,6.35,0,0,0)
    surf.flat(rays,nr=nsil)
    tran.refract(rays,nsil,1.)
    tran.itransform(rays,0,0,0,1.*pi/180,0,0)
    #Go to collimator
    tran.transform(rays,0,0,100,0,0,0)
    surf.flat(rays,nr=1.)
    lenses.collimator6(rays,reverse=True)
    #Go to focus
    tran.transform(rays,0,0,1934.99719-100.,0,0,0)
    surf.flat(rays,nr=1.)
    #Place to AC-508-250
    lenses.AC508_250(rays,reverse=True)
    #Go to WFS location
##    tran.transform(rays,0,0,foc,0,0,0)
##    surf.flat(rays,nr=.1)
    tran.transform(rays,0,0,foc,0,0,0)
    surf.flat(rays,nr=1.)

    #Go to cylindrical field lens
    tran.transform(rays,0,0,-cylz,0,0,0)
    surf.flat(rays,nr=1.)
    tran.transform(rays,0,0,0,0,0,pi/2)
    lenses.LJ1516_L2(rays,reverse=False)
    tran.itransform(rays,0,0,0,0,0,pi/2)
    tran.itransform(rays,0,0,-cylz,0,0,0)
    #Back to WFS
    surf.flat(rays,nr=1.)
    
    return anal.rmsY(rays)

def perfectCyl(rays,align=np.zeros(6)):
    """
    Trace rays from perfect cylinder with potential misalignment
    Assume rays are traced to tangent plane of nominal optic position
    +z points back toward CGH
    Leave with reference frame at tangent plane of nominal surface
    """
    #Apply misalignment
    tran.transform(rays,*align)
    #Trace cylinder
    tran.transform(rays,0,0,220.,0,0,0)
    #Get cylindrical axis in +x direction
    tran.transform(rays,0,0,0,0,0,pi/2)
    surf.cyl(rays,220.,nr=1.)
    tran.reflect(rays)
    tran.itransform(rays,0,0,0,0,0,pi/2)
    tran.itransform(rays,0,0,220.,0,0,0)
    #Go back to nominal tangent plane
    tran.itransform(rays,*align)
    surf.flat(rays,nr=1.)
    
    return

def traceToTestOptic(N):
    """Trace a set of rays from the point source to the nominal
    test optic location
    Return the rays at the plane tangent to the nominal source position.
    """
    #Set up source
    rays = sources.pointsource(.038996,N)
    #Trace through collimator
    tran.transform(rays,0,0,1935.033,0,0,0)
    surf.flat(rays,nr=1.)
    lenses.collimator6(rays)
    #Trace to CGH
    tran.transform(rays,0,0,100.,0,0,0)
    #Apply proper CGH misalignment
    tran.transform(rays,0,0,0,-1.*pi/180,0,0)
    #Trace through CGH
    surf.flat(rays,nr=1.)
    tran.refract(rays,1.,nsil)
    tran.transform(rays,0,0,6.35,0,0,0)
    surf.flat(rays,nr=nsil)
    tran.refract(rays,nsil,1.)
    surf.zernphase(rays,cghcoeff,80.,632.82e-6)
    #Reverse CGH misalignment
    tran.itransform(rays,0,0,0,-1.*pi/180,0,0)
    #Go to line focus
    tran.transform(rays,0,0,0,0,1.*pi/180,0)
    surf.flat(rays,nr=1.)
    tran.transform(rays,0,0,line,0,0,0)
    surf.flat(rays,nr=1.)
    #Go to test optic
    tran.transform(rays,0,0,220.,0,0,0)
    surf.flat(rays,nr=1.)
    #Rotate reference frame so rays impinge toward -z
    tran.transform(rays,0,0,0,0,pi,0)
    
    return rays

def testCollimator():
    """Find out the true focal length of the collimator for
    positioning of the source in raytrace.
    """
    #Create circular plane wave
    rays = sources.circularbeam(25.4*3,100000)
    #Trace through collimator
    lenses.collimator6(rays,reverse=True)
    #Find focus
    foc = anal.analyticImagePlane(rays)
    #Trace to focus and return ray spread and focal distance
    tran.transform(rays,0,0,foc,0,0,0)
    surf.flat(rays,nr=1.)
    return rays
