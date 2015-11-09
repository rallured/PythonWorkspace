import traces.PyTrace as PT
from numpy import *
from matplotlib.pyplot import *
import pdb,time,reconstruct
##from utilities.plotting import mycontour,nanmean
import zernikemod as zmod
import traces.conicsolve as conicsolve

#Import RNR283 figure Zernike coefficients
rnr283coeff = genfromtxt('/home/rallured/Dropbox/AXRO/Alignment/WFS/Simulation/'
                         '140217RNR283Zernikes.txt')
rnr283coeff = rnr283coeff[1:]/1000.
#Import influence Zernike coefficients
infcoeff = genfromtxt('/home/rallured/Dropbox/AXRO/Alignment/WFS/Simulation/'
                      '140217Influence18Zernikes.txt')
infcoeff = infcoeff[1:]/1000.

#Import RNR283 + influence Zernikes
rnrinfcoeff = genfromtxt('/home/rallured/Dropbox/AXRO/Alignment/WFS/Simulation/'
                         '140217Influence18RNR.txt')
rnrinfcoeff = rnrinfcoeff[1:173]/1000.
##rorderri,aorderri = zmod.zmodes(size(rnrinfcoeff))
###Remove coefficients < 1 nm
##ind = where(abs(rnrinfcoeff)>1.e-6)
##rnrinfcoeff = rnrinfcoeff[ind]
##rorderri = rorderri[ind]
##aorderri = aorderri[ind]

#Traces rays from nominal test optic location through the
#system to the image plane
#Misalignments included via arguments, traced using transform
#function
#fieldx = misalignment of field lens in x, fieldy and fieldz analogous
def tracefromtest(fieldx=0.,fieldy=0.,fieldz=0.,imgx=0.,imgy=0.,imgz=0.):
    #Trace to collimator lens
    PT.transform(0.,0.,100.,0.,0.,0.)
    PT.flat() #Get to center of first surface
    PT.lens(2704.01,779.37,10.92,64.01*2,1.64363) #First lens
    PT.transform(0,0,.1,0,0,0) #Gap in achromat
    PT.flat() #Get to center of first surface
    PT.lens(780.87,-1131.72,15.42,64.01*2,1.51501) #Second lens

    #Trace to beamsplitter
    PT.transform(0.,0.,1780.2,0,0,0)
    PT.flat()
    PT.refract(1.,1.51501)#Refract into beamsplitter
    PT.transform(0,0,50,0,0,0)
    PT.flat()
    PT.refract(1.51501,1.)#Refract out of beamsplitter

    #Trace to field lens
    PT.transform(0,0,88.864,0,0,0)
    PT.transform(fieldx,fieldy,fieldz,0,0,0) #Field lens misalignment
    PT.flat()
    PT.lens(205.86,-205.86,5.05,50.,1.51501)
    PT.transform(-fieldx,-fieldy,-fieldz,0,0,0) #Reverse lens misalignment

    #Trace to image plane
    PT.transform(0.,0.,220.085-2.67676768,0,0,0)
    PT.transform(imgx,imgy,imgz,0,0,0) #Include WFS misalignment
    PT.flat()

def fullbeam(num,coeffs,perfect=False,tiltx=0.,tilty=0.,tiltz=0.,**kwgs):
    PT.circularbeam(50.,num) #Define circular beam
    PT.transform(0,0,0,pi,0,0) #Flip beam direction to -z
    PT.transform(0,0,-100,0,0,0) #Move to test optic location

    rorder,aorder = zmod.zmodes(size(coeffs))
    rorder = array(rorder)
    aorder = array(aorder)
##    #Remove coefficients < 1 nm
##    ind = where(abs(coeffs)>1.e-6)
##    coeffs = coeffs[ind]
##    rorder = rorder[ind]
##    aorder = aorder[ind]

    #Trace to RNR283 Zernike
    PT.transform(0,0,0,tiltx,tilty,tiltz)
    if perfect==False:
        PT.zernsurf(coeffs,50.,rorder=rorder,aorder=aorder)
    else:
        PT.flat()
    PT.reflect()
    PT.itransform(0,0,0,tiltx,tilty,tiltz)

    #Trace through rest of system
    tracefromtest(**kwgs)

def centerpointsource(num,**kwgs):
    PT.pointsource(pi/180,num) #Point source at x,y,z=0, points +z

    #Trace through rest of system
    tracefromtest(**kwgs)

#Trace chief and rays on edge of optic
#Rays point to -z and are at +1z
#Purpose is to establish center and radius of
#image for Zernike fitting
#May not be necessary...radius found by maximum without NaNs...
#Use average radius to scale image to some constant radius
#This way, Zernike fits are not really needed...directly compare
#influence functions to one another
def chiefandmarginalrays(**kwgs):
    r = 49.999 #50 mm optic radius
    theta = arange(0,2*pi,pi/4)
    PT.x = r*cos(theta)
    PT.x = append(PT.x,1.e-5)
    PT.y = r*sin(theta)
    PT.y = append(PT.y,1.e-5)
    PT.z = repeat(1.,size(PT.y))
    PT.l = repeat(0.,size(PT.y))
    PT.m = copy(PT.l)
    PT.n = repeat(-1.,size(PT.y))
    PT.ux = copy(PT.l)
    PT.uy = copy(PT.l)
    PT.uz = copy(PT.l)

    PT.zernsurf(rnr283coeff,50.)
    PT.reflect()

    tracefromtest(**kwgs)

def computeinfluence(num,perfect=False,block=False,tiltx=0.,tilty=0.,\
                     tiltz=0.,**kwgs):
    tstart = time.time()
    #Figure out nominal radius
    chiefandmarginalrays()
    nomrad = mean(sqrt(PT.x[:-1]**2+PT.y[:-1]**2))
    
    #Trace chief ray to figure out center of image
    chiefandmarginalrays(**kwgs)
    cx = PT.x[-1]
    cy = PT.y[-1]
    rad = mean(sqrt((PT.x[:-1]-cx)**2+(PT.y[:-1]-cy)**2))

    #Perform actual raytrace without actuator influence
    fullbeam(num,rnr283coeff,perfect=perfect,tiltx=tiltx,tilty=tilty,\
             tiltz=tiltz,**kwgs)
    
    #Subtract off center image coordinate
    PT.x = PT.x - cx
    PT.y = PT.y - cy
    #Rescale to proper radius
    PT.x = PT.x * nomrad/rad
    PT.y = PT.y * nomrad/rad
    
    xang,yang,phase = reconstruct.southwellbin(PT.x,PT.y,PT.l,PT.m,.114,130,130)

    #Perform actual raytrace with actuator influence
    fullbeam(num,rnrinfcoeff,perfect=perfect,**kwgs)
    
    #Subtract off center image coordinate
    PT.x = PT.x - cx
    PT.y = PT.y - cy
    #Rescale to proper radius
    PT.x = PT.x * nomrad/rad
    PT.y = PT.y * nomrad/rad
    xang2,yang2,phase2 = reconstruct.southwellbin(PT.x,PT.y,PT.l,PT.m,.114,130,130)

    #Construct phase and angle arrays for reconstruction of influence
    phaseinf = copy(phase)
    phaseinf[:,:] = 0.
    ind = logical_or(phase==100,phase2==100)
    phaseinf[ind] = 100
    xanginf = copy(xang)
    xanginf = xang2-xang
    yanginf = copy(yang)
    yanginf = yang2-yang
    pdb.set_trace()
    xanginf[ind] = 100
    yanginf[ind] = 100

    #Subtract average tilts
    xanginf[invert(ind)] = xanginf[invert(ind)] - nanmean(xanginf[invert(ind)])
    yanginf[invert(ind)] = yanginf[invert(ind)] - nanmean(yanginf[invert(ind)])

    #Remove portions of measurement to determine effect of alignment marks
    if block is True:
        xanginf[:34,:34] = 100.
        xanginf[:34,98:] = 100.
        xanginf[98:,:34] = 100.
        xanginf[98:,98:] = 100.
        yanginf[:34,:34] = 100.
        yanginf[:34,98:] = 100.
        yanginf[98:,:34] = 100.
        yanginf[98:,98:] = 100.
        phaseinf[:34,:34] = 100.
        phaseinf[:34,98:] = 100.
        phaseinf[98:,:34] = 100.
        phaseinf[98:,98:] = 100.

    #Reconstruct influence wavefront
    influence = reconstruct.reconstruct(xanginf,yanginf,1.e-12,.114,phaseinf)

    #Make invalid pixels NaNs
    ind = where(influence==100)
    influence[ind] = NaN

    print time.time()-tstart

    return influence

#Run alignment scans
def alignscans():
    inom = computeinfluence(10**7)

    #Scan through lateral field lens position
    axialrms = []
    for off in linspace(-5.,5.,10):
        infl = computeinfluence(10**7,fieldx=off)
        axialrms.append(sqrt(nanmean((infl-inom)**2)))
        print axialrms[-1]

    return axialrms

#Put Wolter in full beam system
def fullaperture(num,coeffs,perfect=False,**kwgs):
    PT.circularbeam(50.,num) #Define circular beam
    PT.transform(0,0,0,pi,0,0) #Flip beam direction to -z
    PT.transform(0,0,-100,0,0,0) #Move to test optic location

    #Rotate to proper pitch angle, then translate to focus,
    #place primary, and reverse transformations
    r0 = 220.
    z0 = 8400.
    alpha,p,d,e = conicsolve.woltparam(r0,z0)
    PT.transform(0,0,0,pi/2-alpha+pi,0,0) #Send rays in normal to mirror
    #Wolter misalignment
    PT.transform(0,0,0,pitch,yaw,roll)
    #Go to Wolter focus minus half of mirror and gap
    PT.transform(0,conicsolve.primrad(8476.,r0,z0),-8476.,0,0,0)
    PT.wolterprimary(r0,z0)
    PT.reflect()
    #Go back to center of mirror and reverse misalignment
    PT.transform(0,-conicsolve.primrad(8476.,r0,z0),8476.,-pitch,-yaw,-roll)
    #Rotate to nominal optical axis
    PT.transform(0,0,0,alpha-pi/2,0,0)

    #Trace through rest of system
    tracefromtest(**kwgs)

#Compare theoretical influence (from Zernike fit) to measured influence
def comparetheory(nom=None):
    if nom is None:
        nom = computeinfluence(10**7)
    nom = zmod.stripnans(nom)
    nom = zmod.padimage(nom)/2.

    cx,cy,rad = zmod.locateimage(nom)
    #Loop through and find best center to minimize RMS variation
    theorysurf = zmod.zernsurf(arange(shape(nom)[0]),arange(shape(nom)[1]),\
                               cx,cy,rad,infcoeff)*1000.
    
