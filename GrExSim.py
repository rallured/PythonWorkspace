import PyTrace as PT
from numpy import *
from matplotlib.pyplot import *
import pdb,time,reconstruct
from plotting import mycontour
from plotting import nanmean
import zernikemod as zmod
import conicsolve
import sys
import pickle

#Import efficiencies summed over order -1 through -3
def efficiencyTM(o=None):
##    incidence = arange(.1,3.5+1./60,1./60)
##    i = argmin(abs(incidence-angle))
    f = open('1405052AngScanTMIR.txt','rb')
    d = pickle.load(f)
    f.close()
    wave = d[0]
    eff = d[2]
    order = d[1]
##    wave = wave[i]
##    eff = eff[i]
##    order = order[i]
    effs = zeros(size(unique(wave)))
    l = 0
    for w in unique(wave):
        if o is None:
            effs[l] = sum(eff[where(logical_and(wave==w,\
                            logical_and(order<0,order>-4)))])
        else:
            effs[l] = eff[where(logical_and(wave==w,\
                                            order==o))[0]]
        l = l+1

    return unique(wave),effs

#Same thing for TE
def efficiencyTE(o=None):
##    incidence = arange(.1,3.5+1./60,1./60)
##    i = argmin(abs(incidence-angle))
    f = open('1405052AngScanTEIR.txt','rb')
    d = pickle.load(f)
    f.close()
    wave = d[0]
    eff = d[2]
    order = d[1]
##    wave = wave[i]
##    eff = eff[i]
##    order = order[i]

    effs = zeros(size(unique(wave)))
    l = 0
    for w in unique(wave):
        if o is None:
            effs[l] = sum(eff[where(logical_and(wave==w,\
                            logical_and(order<0,order>-4)))])
        else:
            effs[l] = eff[where(logical_and(wave==w,\
                                            order==o))[0]]
        l = l+1

    return unique(wave),effs


#Analytical focal length determination
def focallengthoffset(rad,xf,zf):
    actual = sqrt((rad-xf)**2+(4750.-zf)**2)
    return 4750.-sqrt(actual**2-rad**2)

#Trace inner grating at r=100mm
#Return location of 0 order focus wrt telescope focus
#This can be used in outer grating trace
def innertrace():
    PT.convergingbeam2(4850.,100.5,103.5,-50.,50.,100000,5.)
    PT.transform(100.,0,4750.,0,pi/2+arctan(100./4750.)-1.5*pi/180,0)
    PT.flat()
    ind = where(logical_and(abs(PT.x)<50.,abs(PT.y)<50.))
    PT.vignette(ind=ind)

    PT.reflect()
    PT.transform(0,0,0,0,-pi/2-arctan(100./4750.)+1.5*pi/180,0)
    PT.transform(-100.,0,-4750.,0,0,0)
    PT.flat()

    zoffset = PT.findimageplane2(3.,500.)
    PT.transform(0,0,zoffset,0,0,0)
    cx = mean(PT.x)
    cy = mean(PT.y)
    PT.transform(cx,cy,0,0,0,0)

    return cx,cy,zoffset

##xf,yf,zf = innertrace()
##focaloffset = -27.9797979797

#Trace outer grating
#Optimize focal length of converging beam and pitch
#of grating to produce 0 order matched to inner grating
def traceouter(rad,focaloffset,xf=None,yf=None,zf=None):
    #Set converging beam to outer grating
    PT.convergingbeam2(4850.+focaloffset,662.,665.5,-50.,50.,10000000,5.)
    if xf is None:
        xf,yf,zf = innertrace()
    #Create vector pointing to telescope focus
    xt = -650.
    zt = -4750.-focaloffset
    mag = sqrt(xt**2+zt**2)
    xt = xt/mag
    zt = zt/mag
    #Create vector pointing to inner 0 order
    x0 = xf - 650.
    z0 = zf - 4750.
    mag = sqrt(x0**2+z0**2)
    x0 = x0/mag
    z0 = z0/mag
    #Dot them to determine angular deviation
    ang = arccos(x0*xt + z0*zt)
    #Determine required pitch angle to deviate beam to 0 order focus
    PT.transform(650.,0,4750.+focaloffset,0,pi/2+arctan(650./\
                                (4750.+focaloffset))-ang/2,0)
    PT.flat()
    ind = where(logical_and(abs(PT.x)<50.,abs(PT.y)<50.))
    PT.vignette(ind=ind)

    #Reflect and trace to inner grating focal plane
    #Return HEW and required pitch
    PT.reflect()
    PT.transform(0,0,0,0,-pi/2-arctan(650./\
                            (4750.+focaloffset))+ang/2,0)
    PT.transform(-650.+xf,0,-4750.+zf,0,0,0)
    PT.flat()
    rho = sqrt(PT.x**2+PT.y**2)
    return median(rho),ang/2

#Scan through focal offsets and determine best
#alignment for outer grating
def alignouter(xf,yf,zf):
    offsets = linspace(-30.,0.,31)
    hew = []
    ang = []
    for o in offsets:
        h,a=traceouter(o,xf=xf,yf=yf,zf=zf)
        hew.append(h)
        ang.append(a)
        sys.stdout.write(str(o)+'\r')
        sys.stdout.flush()
    #Determine best offset
    clf()
    plot(offsets,hew)
    return offsets[argmin(hew)],ang[argmin(hew)]

###Figure out z(r)
##def offsetdependence():
##    rad = 

#Trace inner
#Put correct yaw for blaze angle, hub for throw,
#Set nominal 160 nm d spacing at center of grating
def simulateinner(xf,yf,zf,order,wave,rotate=0.):
    #Set up inner grating converging beam
    PT.convergingbeam2(4850.,100.5,103.5,-50.,50.,1000,5.)
    PT.transform(100.,0,4750.,0,pi/2+arctan(100./4750.)-1.5*pi/180,0)
    PT.flat()
    ind = where(logical_and(abs(PT.x)<50.,abs(PT.y)<50.))
    PT.vignette(ind=ind)

    #Compute throw
    throw = sqrt((100.-xf)**2+(4750.-zf)**2)
    inc = 1.5*pi/180
    hubdist = throw/cos(inc)
    #Compute yaw
    blaze = 18.*pi/180
    phi0 = arccos(sin(inc)/cos(blaze))
    yaw = pi/2-arctan(sin(blaze)*(1./tan(phi0)))
    yaweff = -yaw+pi/2
    #Compute dispersion
    dpermm = 160./hubdist
    alpha = sin(blaze)*cos(phi0)+3.8/160.*sin(yaweff)
    beta = sin(phi0)-3.8/160.*cos(yaweff)
    throw = throw*cos(inc)
    disp = (throw/beta)*(sin(yaw)/160.)+\
           (-throw*alpha/beta**2)*(-cos(yaw)/160.)

    #Rotate to point y axis in groove direction
##    yaweff = 0.
    PT.transform(0,0,0,0,0,-pi/2-yaweff)
    #Diffract and trace to focal plane
    PT.radgrat(hubdist,dpermm,order,wave)
    #Rotate back to original coordinate system
    PT.transform(0,0,0,0,0,pi/2+yaweff)
    PT.transform(0,0,0,0,-pi/2-arctan(100./4750.)+inc,0)
    #Go to focal plane
    PT.transform(-100+xf,0,-4750.+zf,0,rotate,0)
    PT.flat()
    print 'X: ' + str(nanmean(PT.x))
    print 'Y: ' + str(nanmean(PT.y))

    #Return ray intercepts in dispersion direction (y)
    plot(PT.y,PT.x,'.')
    return PT.y

#Trace outer grating
#Put correct yaw for blaze angle, hub for throw,
#and dispersion for outer grating
def simulateouter(xf,yf,zf,focaloffset,order,wave,rotate=0.,doff=0.):
    #Set converging beam to outer grating
    PT.convergingbeam2(4850.+focaloffset,662.,665.5,-172.,172.,100000,5.)
##    PT.convergingbeam2(4850.+focaloffset,662.,665.5,-50.,50.,1000,5.)
    if xf is None:
        xf,yf,zf = innertrace()
    #Create vector pointing to telescope focus
    xt = -650.
    zt = -4750.-focaloffset
    mag = sqrt(xt**2+zt**2)
    xt = xt/mag
    zt = zt/mag
    #Create vector pointing to inner 0 order
    x0 = xf - 650.
    z0 = zf - 4750.
    mag = sqrt(x0**2+z0**2)
    x0 = x0/mag
    z0 = z0/mag
    #Dot them to determine angular deviation
    inc = arccos(x0*xt + z0*zt)/2
    #Determine required pitch angle to deviate beam to 0 order focus
    PT.transform(650.,0,4750.+focaloffset,0,pi/2+arctan(650./\
                                (4750.+focaloffset))-inc,0)
    PT.flat()
    ind = where(logical_and(abs(PT.x)<50.,abs(PT.y)<170.))
    PT.vignette(ind=ind)
    
    #Compute throw
    throw = sqrt((650.-xf)**2+(4750.-zf)**2)
    hubdist = throw/cos(inc)
    #Compute yaw
    blaze = 18.*pi/180
    phi0 = arccos(sin(inc)/cos(blaze))
    yaw = pi/2-arctan(sin(blaze)*(1./tan(phi0)))
    #Compute dispersion
    disp = 29.69824541 #from inner trace
    alpha = sin(blaze)*cos(phi0)
    beta = sin(phi0)
    yaweff = -yaw+pi/2
    throwproj = throw*cos(inc)
    d = throwproj*((1/beta)*(sin(yaw)/disp)+\
                   (alpha/beta**2)*(cos(yaw)/disp))+doff
    #Required d to put 3.8 nm beam at -112.8501
##    goal = -148.497411416e-3
##    d = 5.0*((alpha*throwproj+beta*goal)*cos(yaw)+\
##                beta*throwproj*sin(yaw))/(beta**2*goal)
    dpermm = d/hubdist

##    dpermm = 160./hubdist

    #Rotate to point y axis in groove direction
    PT.transform(0,0,0,0,0,-pi/2-yaweff)
    #Diffract and trace to focal plane
    PT.radgrat(hubdist,dpermm,order,wave)
    #Rotate back to original coordinate system
    PT.transform(0,0,0,0,0,pi/2+yaweff)
    PT.transform(0,0,0,0,-pi/2-arctan(650./\
                            (4750.+focaloffset))+inc,0)
    #Go to focal plane
    PT.transform(-650+xf,0,-4750.+zf,0,rotate,0)
    PT.flat()

    #Return ray intercepts in dispersion direction (y)
    plot(PT.y,PT.x,'.')
    return PT.y

#Determine resolution for combined gratings
def combinedres(xf,yf,zf,focaloffset,order,wave,rotate=0.,doff=0.):
    inner = simulateinner(xf,yf,zf,order,wave,rotate=rotate)
    inner = inner[invert(isnan(inner))]
    print 'Inner Mean: ' + str(nanmean(inner))
    outer = simulateouter(xf,yf,zf,focaloffset,order,wave,\
                          rotate=rotate,doff=doff)
    outer = outer[invert(isnan(outer))]
    print 'Outer Mean: ' + str(nanmean(outer))
    #Take equal number of counts
    sz = min((size(inner),size(outer)))
    combined = concatenate((inner[:sz],outer[:sz]))
    combined = combined[invert(isnan(combined))]
    width = median(abs(combined-nanmean(combined)))*2
    width = sqrt(width**2+(1./60**2*pi/180*5000.)**2)
    res = nanmean(combined)/width

    print 'Width: ' + str(width)

    return res

#Plot resolution vs. wavelength for order 1
def resolution1(xf,yf,zf,focaloffset):
    wave = linspace(2.5,5.0,30)
    res = []
    for w in wave:
        res.append(combinedres(xf,yf,zf,focaloffset,1.,w,\
                               rotate=7.*pi/180,doff=.1))
    return array(wave),array(res)

#Plot resolution vs. wavelength for order 2
def resolution2(xf,yf,zf,focaloffset):
    wave = linspace(.5,3.0,30)
    res = []
    for w in wave:
        res.append(combinedres(xf,yf,zf,focaloffset,2.,w,\
                               rotate=7.*pi/180,doff=.1))
    return array(wave),array(res)

#Plot resolution vs. wavelength for order 1
def resolution3(xf,yf,zf,focaloffset):
    wave = linspace(.5,3.0,30)
    res = []
    for w in wave:
        res.append(combinedres(xf,yf,zf,focaloffset,3.,w,\
                               rotate=7.*pi/180,doff=.1))
    return array(wave),array(res)

#Scan focal plane offset and rotation angle
def scanfocalplane(xf,yf,zf,focaloffset):
    #Set up offset range and rotation angle range
    rotation = linspace(4.,7.,20)*pi/180.
    res25 = zeros(size(rotation))
    res50 = copy(res25)
    for i in range(size(rotation)):
        c,res = combinedres(xf,yf,zf,focaloffset,1.,2.5,rotate=rotation[i],doff=.1)
        res25[i] = res
        c,res = combinedres(xf,yf,zf,focaloffset,1.,4.3,rotate=rotation[i],doff=.1)
        res50[i] = res
        print rotation[i]
        sys.stdout.flush()
    pdb.set_trace()

#Test radial gratings
def testrad(wave):  
    #Set up input beam
    PT.x = zeros(10)
    PT.y = zeros(10)
    PT.z = zeros(10)
    PT.l = zeros(10)
    PT.m = repeat(cos(1.5*pi/180),10)
    PT.n = repeat(sin(1.5*pi/180),10)
    PT.ux = zeros(10)
    PT.uy = zeros(10)
    PT.uz = zeros(10)

    #Diffract
    PT.radgrat(8000.,160./8000.,1,wave)
    pdb.set_trace()
    #Rotate
    PT.transform(0,0,0,-pi/2,0,0)
    PT.transform(0,0,8000.,0,0,0)
    PT.flat()

    plot(PT.x,PT.y,'.')

#Effect of subaperture angle
def subaperturetest(angrange):
    #Radius vector
    r = random.uniform(low=-.602/2,high=.602/2,size=1000)

    #Loop through subaperture angle
    deviation = []
    for ang in angrange:
        theta = random.uniform(low=-ang/2.,high=ang/2.,size=1000)
        x = r*sin(theta*pi/180)
        deviation.append(std(x))
    return array(deviation)

#Optimal SPO height calculation
#1: For given R and phi, compute path difference between
#focus and 1st order (have to match diffraction at middle of arc)
#2: Trace from grating to SPO join plane then to offset position.
#Determine height difference and new R'
#3: Convert to cartesian coordinates x', y',
#associate height and nominal R to grid of x',y'
#4: Interpolate to regular grid of x',y', make contour plots
def spoheights():
    #Define nominal focus and 0 order position
    #0 Order from Zemax sensitivity simulation, optimized for 1.5 deg
    #incidence at R=100mm
    x0 = 248.816477631
    z0 = .88376754

    #Loop through R and phi
    rrange = linspace(100.,660.,1000)
    phirange = linspace(-15.,15.,1000)
    h = zeros((size(rrange),size(phirange)))
    ract = copy(h)
    rnom = copy(h)
    for ri in range(size(rrange)):
        for pi in range(size(phirange)):
            r = rrange[ri]
            phi = phirange[pi]
            x = r*cos(phi*pi/180)
            y = r*sin(phi*pi/180)
            phi = phirange[pi]
            #Compute path length difference
            focdist = sqrt(r**2+4750.**2)
            zerodist = sqrt((x-x0)**2+y**2+(4750-z0)**2)
            pathdiff = focdist - zerodist

            #Trace from grating up to required SPO position
            coneangle = arctan(r/4750.)
            rnom[ri,pi] = r + 150.*tan(coneangle)
            ract[ri,pi] = rnom[ri,pi] + sin(coneangle)*pathdiff            
            h[ri,pi] = cos(coneangle)*pathdiff
    pdb.set_trace()
