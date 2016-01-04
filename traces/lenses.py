#Module to define singlet and doublet lenses
#Explicitly put in commonly used lenses such as Thorlabs doublets
import surfaces as surf
import transformations as tran

def singlet(rays,r1,r2,thick,nl,reverse=False):
    """Trace a spherical singlet lens. Assume reference frame
    is +z toward optical axis, xy plane tangent to first surface.
    Positive R indicates convex surface for both surfaces.
    """
    if reverse is True:
        r1,r2 = r2,r1
    #Trace to first surface
    tran.transform(rays,0,0,r1,0,0,0)
    surf.sphere(rays,r1,nr=1.)
    #Refract into material
    tran.refract(rays,1.,nl)
    #Trace to second surface
    tran.transform(rays,0,0,-r1+thick-r2,0,0,0)
    surf.sphere(rays,r2,nr=nl)
    tran.transform(rays,0,0,r2,0,0,0)
    #Refract out of surface
    tran.refract(rays,nl,1.)
    return
    

def lens(rays,r1,r2,thick,d,nl,reverse=False):
    """Trace lens, first surface center is coincident with xy plane
    thickness extends in positive z direction
    Rays should be traced to this plane before calling lens
    to ensure rays hit correct spherical intersection
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if reverse is True:
        r1,r2 = -r2,-r1
    #Trace to first surface
    if r1 != 0:
        transform(rays,0.,0.,r1,0.,0.,0.) #Go to first
                                     #center of curvature, r1>0->convex
        sphere(rays,abs(r1)) #Trace to first spherical surface
    rho = np.sqrt(x**2 + y**2)/(d/2)
    ind = np.where(rho<=1.) #Remove rays with radius beyond that of lens
    vignette(rays,ind=ind)
    refract(rays,1.,nl) #Refract into lens
    
    transform(rays,0.,0.,-r1+thick,0.,0.,0.) #Go to center of second surface
    flat(rays)
    if r2 != 0:
        transform(rays,0.,0.,r2,0,0,0) #Go to center of curvature
        sphere(rays,abs(r2)) #Trace to second spherical surface
    rho = np.sqrt(x**2 + y**2)/(d/2)
    ind = np.where(rho<=1.) #Remove rays with radius beyond that of lens
    vignette(rays,ind=ind)
    refract(rays,nl,1.) #Refract out of lens
    transform(rays,0.,0.,-r2,0,0,0) #Transform to xy plane tangent to
                              #center of second surface

    #Rays are left at second surface, pointing out of lens in correct direction
    return

def cyllens(rays,r1,r2,thick,width,height,nl,reverse=False):
    """Cylindrical lens, same principle as with standard lens
    """
    opd,x,y,z,l,m,n,ux,uy,uz
    if reverse is True:
        r1,r2 = -r2,-r1
    #Trace to first surface
    if r1 != 0:
        transform(rays,0.,0,r1,0.,0.,0.) #Go to first
                                     #center of curvature, r1>0->convex
        cyl(rays,abs(r1)) #Trace to first spherical surface
    ind = logical_and(x<width/2,y<height/2) #Remove rays outside of cylinder
    vignette(rays,ind=ind)
    refract(rays,1.,nl) #Refract into lens

    
    transform(rays,0.,0,-r1+thick,0.,0.,0.) #Go to center of second surface
    flat(rays)
    if r2 != 0:
        transform(rays,0.,0,r2,0,0,0) #Go to center of curvature
        cyl(rays,abs(r2)) #Trace to second spherical surface
    ind = logical_and(x<width/2,y<height/2) #Remove rays outside of cylinder
    vignette(rays,ind=ind)
    refract(rays,nl,1.) #Refract out of lens
    transform(rays,0.,0,-r2,0,0,0) #Transform to xy plane tangent to
                              #center of second surface

    #Rays are left at second surface, pointing out of lens in correct direction
    return


####### Lenses ##########
def collimator6(rays,reverse=False):
    """
    Traces through the six inch collimator from Cumberland.
    R1=1124.
    R2=9324.
    Standard orientation is collimation of point source
    Reverse is focusing of plane wave
    """
    singlet(rays,9324.,1124.,20.,1.5150885,reverse=reverse)
    return

def LJ1653L2(rays,reverse=False):
    cyllens(rays,103.36,0,4.09,30.,60.,1.51501,reverse=reverse)
    return

def LJ1629L2(rays,reverse=False):
    cyllens(rays,77.52,0,4.46,30.,60.,1.51501,reverse=reverse)
    return

def AC254_400_A(rays,reverse=False):
    if reverse is False:
        lens(rays,738.5,181.55,2.,12.7*2,1.64363)
        lens(rays,181.55,-219.8,4.,12.7*2,1.51501)
    else:
        lens(rays,181.55,-219.8,4.,12.7*2,1.51501,reverse=True)
        lens(rays,738.5,181.55,2.,12.7*2,1.64363,reverse=True)
    return

def AC508_200_A(rays,reverse=False):
    if reverse is False:
        lens(rays,109.86,-93.110,8.5,50.8,1.51501)
        lens(rays,-93.110,-376.25,2.,50.8,1.64363)
    else:
        lens(rays,-93.110,-376.25,2.,50.8,1.64363,reverse=True)
        lens(rays,109.86,-93.110,8.5,50.8,1.51501,reverse=True)
    return

def cylNull(rays,reverse=False):
    if reverse is False:
        transform(rays,0,0,0,np.pi/2,0,0)
        cylconic(rays,.007626,-.575)
        refract(rays,1.,1.51501)
        transform(rays,0,0,0,-np.pi/2,0,0)
        transform(rays,0,0,50,0,0,0)
        flat(rays,)
        refract(rays,1.51501,1.)
    else:
        refract(rays,1.,1.51501)
        transform(rays,0,0,50,0,0,0)
        transform(rays,0,0,0,np.pi/2,0,0)
        cylconic(rays,-.007626,-.575)
        refract(rays,1.51501,1.)
        transform(rays,0,0,0,-np.pi/2,0,0)
    return
