import traces.PyTrace as PT
import numpy as np
import traces.conicsolve as con
import pdb,sys
from traces.axro.WSverify import traceChaseParam
import matplotlib.pyplot as plt
import time

#Need to determine resolution and effective area
#as a function of shell radius
#Loop through all ~260 shells and determine vignetting
#and resolution functions vs. pointing error
#Will also need reflectivity vs. theta

#IMD Ir reflectivity at 1 keV
irang,irref = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO'
                            '/WSTracing/'
                            '150504IrReflectivity.txt',comments=';'))
#CXC Ir optical constants
ener,delta,beta,alpha,gamma = np.transpose(np.genfromtxt('/home/rallured'
                        '/Dropbox/AXRO/WSTracing/chandraConstants.txt')[2:])
#Paul's thermal filter transmission


def reflectivityIr(ang):
    """Return reflectivity with 0.5 nm RMS roughness
    calculated with IMD (1 keV)"""
    return irref[np.argmin(np.abs(irang-ang))]

def CXCreflIr(ang,energy,rough):
    """Return reflectivity with 0.5 nm RMS roughness
    calculated using Fresnel equations, Chandra
    optical constants, and "Strehl" factor (Nevot-Croce)
    Energy supplied in eV
    Roughness in RMS nm
    """
    #Get proper optical constants
    ind = np.argmin(abs(ener-energy/1000.))
    b = beta[ind]
    d = delta[ind]
    n = 1 - d + 1j*b
    n2 = abs(n)**2
    #Compute reflectivity in each polarization plane
    #Return mean value
    Rp = abs(n*n*np.sin(ang)-np.sqrt(n*n-np.cos(ang)**2))**2/\
         abs(n*n*np.sin(ang)+np.sqrt(n*n-np.cos(ang)**2))**2
    Rs = abs(np.sin(ang)-np.sqrt(n*n-np.cos(ang)**2))**2/\
         abs(np.sin(ang)+np.sqrt(n*n-np.cos(ang)**2))**2
    R = np.mean([Rp,Rs],axis=0)
    wave = 1240./energy #wavelength in nm
    k = 2*np.pi/wave
    strehl = np.exp(-4*k**2*np.sin(ang)**2*rough**2)

    return R*strehl

def traceWSShell(num,theta,r0,z0,phigh,plow,shigh,slow,\
                 energy,rough,chaseFocus=False,bestFocus=False):
    """Trace a WS mirror pair with 10 m focal length and
    mirror axial cutoffs defined by phigh,plow
    """
    #Define annulus of source rays
    a,p,d,e = con.woltparam(r0,z0)
    r1 = PT.wsPrimRad(plow,1.,r0,z0)#np.tan(a/2.)*(plow-10000.) + r0
    r2 = PT.wsPrimRad(phigh,1.,r0,z0)#np.tan(a/2.)*(phigh-10000.) + r0
    PT.annulus(r1,r2,num)
    PT.transform(0,0,0,np.pi,0,0)
    PT.transform(0,0,z0,0,0,0)

    #Trace to primary
    PT.wsPrimary(r0,z0,1.)
    #Handle vignetting
    PT.vignette()
    ind = np.logical_and(PT.z<phigh,PT.z>plow)
    PT.vignette(ind=ind)
    #Vignette rays hitting backside of mirror
    dot = PT.l*PT.ux+PT.m*PT.uy+PT.n*PT.uz
    ind = dot < 0.
    PT.vignette(ind=ind)
    #If all rays are vignetted, return
    if np.size(PT.x) < 1:
        return 0.,0.,0.
    #Apply pointing error
    PT.l = PT.l + np.sin(theta)
    PT.n = -np.sqrt(1 - PT.l**2)
    #Reflect
    PT.reflect()
    #Compute mean incidence angle for reflectivity
    ang = np.abs(np.mean(np.arcsin(dot))) #radians
    refl1 = CXCreflIr(ang,energy,rough)
    #Total rays entering primary aperture
    N1 = np.size(PT.x)

    #Trace to secondary
    PT.wsSecondary(r0,z0,1.)
    #Vignette anything that did not converge
    PT.vignette()
    #Vignette anything outside the physical range of the mirror
    ind = np.logical_and(PT.z>slow,PT.z<shigh)
    PT.vignette(ind=ind)
    #Vignette anything hitting the backside
    dot = PT.l*PT.ux+PT.m*PT.uy+PT.n*PT.uz
    ind = dot < 0.
    PT.vignette(ind=ind)
    if np.size(PT.x) < 1:
        return 0.,0.,0.
    PT.reflect()
    #Compute mean incidence angle for reflectivity
    ang = np.abs(np.mean(np.arcsin(dot))) #radians
    refl2 = CXCreflIr(ang,energy,rough)

    #Trace to focal plane
    PT.flat()

    #Find Chase focus
    delta = 0.
    if chaseFocus or bestFocus:
        delta = .0625*(psi+1)*(r**2*L1/z0**2)*(1/np.tan(alpha))**2
        PT.transform(0,0,delta,0,0,0)
        PT.flat()

    #Find best focus
    delta2 = 0.
    delta3 = 0.
    if bestFocus:
        try:
            delta2 = PT.findimageplane(20.,100)
            PT.transform(0,0,delta2,0,0,0)
            delta3 = PT.findimageplane(1.,100)
            PT.transform(0,0,delta3,0,0,0)
        except:
            pdb.set_trace()
        PT.flat()

    return refl1*refl2

def evaluateShell(theta,r0):
    """Compute vignetting factor and HPD as a function of
    pointing error. Supply pointing errors theta and intersection
    radius of shell. Assumption is 200 mm long segments.
    """
    hpd = np.zeros(np.size(theta))
    vign = np.copy(hpd)
    refl = np.copy(hpd)
    a,p,d,e = con.woltparam(r0,10000.)
    for t in theta:
        hpd[t==theta],vign[t==theta],refl[t==theta] = \
            traceWSShell(1000,t,r0,10225.,10025.,9975.,9775.)

    return hpd,vign,refl

def SXperformance(theta,energy,rough,bestsurface=False):
    """Go through a SMART-X prescription file and compute
    area weighted performance for a flat focal plane
    """
    #Load in rx data
##    rx = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO/WSTracing/'
##        'mirror-design-260sh-200mmlong-040mmthick'
##        '-3mdiam-10mfl-10arcmin-fov-planarIntersept032713.csv',\
##                       delimiter=','))
    rx = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO/WSTracing/'
                                    '150528_Pauls_Rx.csv',delimiter=','))
    geo = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO/WSTracing/'
                                     'geometric_transmission_102711.txt'))
    therm = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO/'
                                       'WSTracing/thermal_shield_transmission_102711.txt'))
    f = np.sqrt(rx[1][-1]**2+10000.**2)
    z = np.sqrt(f**2-rx[1]**2)
##    z = np.repeat(10000.,np.size(rx[1]))
##    ind = rx[0] > 210.
##    rx = rx[:,ind]
    Ns = np.shape(rx)[1]
    pdb.set_trace() 

    #Loop through and compute a resolution and a weight for each shell
    hpdTelescope = np.zeros(np.size(theta))
    rmsTelescope = np.zeros(np.size(theta))
    delta = np.zeros(np.size(theta))
    #fig = plt.figure()
    for t in theta[:]:
        xi = np.array([])
        yi = np.array([])
        l = np.array([])
        m = np.array([])
        n = np.array([])
        weights = np.array([])
        #plt.clf()
        tstart = time.time()
        for s in np.arange(0,Ns):
            if geo[1][s] > 0.:
                sys.stdout.write('Shell: %03i \r' % s)
                sys.stdout.flush()
                r = traceWSShell(1000,t,rx[1][s],z[s],z[s]+225.,z[s]+25.,\
                                     z[s]-25.,z[s]-225.,energy,rough)
                r = r*geo[1][s]*rx[9][s] #Reflectivity*area*alignmentbars*vign
                #Account for thermal shield in shells 220-321
                if s > 219:
                    r = r * therm[1][np.abs(energy/1000.-therm[0]).argmin()]
                r = np.repeat(r,np.size(PT.x))
                weights = np.append(weights,r)
                xi = np.append(xi,PT.x)
                yi = np.append(yi,PT.y)
                l = np.append(l,PT.l)
                m = np.append(m,PT.m)
                n = np.append(n,PT.n)
                #plt.plot(PT.x,PT.y,'.')
        print time.time()-tstart

        #Have list of photon positions and weights
        #Need to compute centroid and then FoM
        #Normalize weights
        weights = weights/np.sum(weights)
        PT.x = np.array(xi,order='F')
        PT.y = np.array(yi,order='F')
        PT.z = np.zeros(np.size(xi)).astype('float')
        PT.l = np.array(l,order='F')
        PT.m = np.array(m,order='F')
        PT.n = np.array(n,order='F')
        PT.ux = np.zeros(np.size(xi)).astype('float')
        PT.uy = np.zeros(np.size(xi)).astype('float')
        PT.uz = np.zeros(np.size(xi)).astype('float')
        if bestsurface:
            PT.transform(0,0,.25,0,0,0)
            delta[t==theta] = PT.findimageplane(.5,201,weights=weights)
            PT.transform(0,0,delta[t==theta],0,0,0)
            PT.flat()
        #Compute FoM
        rmsTelescope[t==theta] = PT.rmsCentroid(weights=weights)/10000.
        hpdTelescope[t==theta] = PT.hpd(weights=weights)/10000.
        print hpdTelescope[t==theta],rmsTelescope[t==theta]
       
        pdb.set_trace()

    return hpdTelescope,rmsTelescope,delta

def sphericalNodes(rin,z0,fov,Nshells,N):
    """This function will iteratively scan node positions
    about a sphere around the focus. Node will start in obvious
    vignetting position. Extreme rays will be traced including
    FoV. Node will be nudged outward until vignetting no longer
    occurs. Node will then be moved by the designated mechanical
    gap. Then the next node is traced in the same fashion.
    Assumptions: 50 mm symmetric gap
    """
    #Bookkeeping parameters
    f = np.sqrt(rin**2+z0**2)
    fov = fov/60.*np.pi/180. #fov to radians
    zlist = []
    rlist = []
    for i in range(Nshells):
        #Starting radius for next shell node
        rstart = PT.wsPrimRad(z0+225.,1.,rin,z0)
        #Reduce rstart until vignetting is reached
        flag = 0
        while flag==0:
            zstart = np.sqrt(f**2-rstart**2)
            #Set up rays
            r1 = PT.wsPrimRad(zstart+25.,1.,rstart,zstart)
            r2 = PT.wsPrimRad(zstart+225.,1.,rstart,zstart)
            PT.pointsource(0.,N)
            PT.z = np.repeat(10500.,N)
            PT.x = np.linspace(r1,r2,N)
            PT.n = np.repeat(-1.,N)
            #Perform trace and add FoV deflections to rays
            PT.wsPrimary(rstart,zstart,1.)
            PT.l = np.repeat(np.sin(fov),N)
            PT.n = -np.sqrt(1 - PT.l**2)
            #Verify that rays do not hit prior primary
            PT.wsPrimary(rin,z0,1.)
            if np.sum(PT.z<z0+225.) != 0:
                #Ray has hit
                print 'Ray hits prior primary!'
                flag = 1
            #Verify that rays do not hit prior secondary
            PT.wsPrimary(rstart,zstart,1.)
            PT.reflect()
            PT.wsSecondary(rstart,zstart,1.)
            PT.reflect()
            PT.wsSecondary(rin,z0,1.)
            if np.sum(PT.z > z0-225.) != 0:
                print 'Ray hits prior secondary!'
                flag = 1
            #Look at other deflection
            PT.pointsource(0.,N)
            PT.z = np.repeat(10500.,N)
            PT.x = np.linspace(r1,r2,N)
            PT.n = np.repeat(-1.,N)
            #Perform trace and add FoV deflections to rays
            PT.wsPrimary(rstart,zstart,1.)
            PT.l = np.repeat(-np.sin(fov),N)
            PT.n = -np.sqrt(1 - PT.l**2)
            #Verify that rays do not hit prior primary
            PT.wsPrimary(rin,z0,1.)
            if np.sum(PT.z<z0+225.) != 0:
                #Ray has hit
                print 'Ray hits prior primary!'
                flag = 1
            #Verify that rays do not hit prior secondary
            PT.wsPrimary(rstart,zstart,1.)
            PT.reflect()
            PT.wsSecondary(rstart,zstart,1.)
            PT.reflect()
            PT.wsSecondary(rin,z0,1.)
            if np.sum(PT.z > z0-225.) != 0:
                print 'Ray hits prior secondary!'
                flag = 1
            if flag==0:
                rstart = rstart - .01 #Take off 10 microns
##                sys.stdout.write(str(rstart)+'\n')
##                sys.stdout.flush()
        #Vignetting has been reached, append rstart and zstart
        #to list of node positions
        rlist.append(rstart)
        zlist.append(zstart)
        rin = rstart
        z0 = zstart
    return rlist,zlist


