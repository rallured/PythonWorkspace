import numpy as np
from numpy import sin,cos,exp,sqrt,pi,tan
import matplotlib.pyplot as plt
import pdb,sys
import traces.grating as grat
import utilities.plotting as plotting
from traces.axro.SMARTX import CXCreflIr
from scipy import interpolate
import utilities.imaging.man as man
import utilities.transformations as tr
import utilities.imaging.fitting as fit
import astropy.convolution as conv

import traces.PyTrace as PT

import traces.analyses as anal
import traces.surfaces as surf
import traces.transformations as tran
import traces.sources as sources

#Load CCD QE data and define interpolation function
ccd = np.genfromtxt('/home/rallured/Dropbox/'
                    'Arcus/Raytrace/151002_CCDQE.csv',delimiter=',')
ccd = np.transpose(ccd)
ccdQE = interpolate.interp1d(1239.8/ccd[0],ccd[1],kind='linear')
#Load efficiency data and return interpolation function
w = np.genfromtxt('/home/rallured/Dropbox/'
                'Arcus/Raytrace/Wave.txt')
o = np.genfromtxt('/home/rallured/Dropbox/'
                  'Arcus/Raytrace/Order.txt')
eff = np.genfromtxt('/home/rallured/Dropbox/'
                    'Arcus/Raytrace/Eff.txt')

#Nominal OPG design
bestfoc = -11733.773388505848
#Standard deviation in microns for allowable
#contribution to LSF of everything aside from
#design aberrations (i.e. SPO quality, alignments, everything else)
marg = 0.048027499983957431


def gratEff(order):
    return interpolate.interp1d(w,eff[:,o==order],kind='linear',axis=0)

##How to investigate background?
##Turn wavelength into a vector. Draw from distribution Andy sent you.
##Also apply random wavevector kick. But how big? This should be limited by
##geometric stray light approximations. Then trace things through and
##return list of background photons at the focal plane. This shouldn't
##be too bad.

def meritFunction():
    """Go through the science cases, trace the relevant lines,
    sum the observation times
    """
    #Set up system
    yaw = grat.blazeYaw(1.5*pi/180,2.4,3,160.)
    bestfoc = sec.traceSector(300.,800.,12e3,100,span=375.,\
                              blazeYaw=yaw,bestFocus=None,order=-3,wave=2.4)
    marg = sec.traceSector(300.,800.,12e3,100,span=375.,blazeYaw=yaw,\
                           bestFocus=bestfoc,order=-3,wave=2.4,\
                           findMargin=True)
    
    #Case G1-1, column densities as function of radial distance
    #from cluster center
    wave = 2.16
    o = np.arange(1,8)
    o = o[(4.4<wave*o) & (wave*o<8.8)]

    return None
    


def investigateSector(Rin,Rout,F,N,wave,span=20.,d=.605,t=.775,gap=50.,\
                inc=1.5*pi/180,l=95.,\
                blazeYaw=0.,order=-1,marg=0.):
    #Get axial offset at design wavelength
    bestFoc = traceSector(Rin,Rout,F,N,span=span,d=d,t=t,gap=gap,\
                inc=inc,l=l,order=-3,\
                blazeYaw=blazeYaw,wave=2.4)
    
    #Get axial offset at nominal focus - this sets the
    #focal plane tilt angle
##    bestFoc0 = traceSector(Rin,Rout,F,N,span=span,d=d,t=t,gap=gap,\
##                inc=inc,l=l,order=-1,\
##                blazeYaw=blazeYaw,wave=4.9)[0]
    
    #Get grating efficiency function
    geff = gratEff(order)
    
    #Set up loops to investigate resolution and
    #effective area as a function of wavelength
    res = np.zeros(np.size(wave))
    area = np.copy(res)
    for i in range(np.size(wave)):
        b,r,a = 0.,0.,0.
        if geff(wave[i]) > 0.:
            r,a = traceSector(Rin,Rout,F,N,span=span,d=d,\
                            t=t,gap=gap,inc=inc,l=l,\
                            wave=wave[i],blazeYaw=blazeYaw,\
                            bestFocus=bestFoc,order=order,marg=marg)[:2]
        #Add in grating efficiency and CCD QE
##        a = a * geff(wave[i]) * ccdQE(wave[i])
        res[i] = r
        area[i] = a
        sys.stdout.write('Wave: %.2f, Area: %.2f, Res: %.2f\r'\
                         % (wave[i],area[i],res[i]))
        sys.stdout.flush()
        
    return res,area

def traceSector(Rin,Rout,F,N,span=20.,d=.605,t=.775,gap=50.,\
                inc=1.5*pi/180,l=95.,bestFocus=None,order=0,\
                blazeYaw=0.,wave=1.,marg=0.,findMargin=False):
    """Trace Arcus sector
    """
    #Trace parameters
    R = np.arange(Rin,Rout,t) #Vector of shell radii
    tg = .25*np.arctan((R+d/2)/F) #Primary graze angles
    L = d/tan(tg) #Primary length
    #L = 4*F*d/R #Vector of mirror lengths
    M = np.size(R) #Number of shells

    #Create focal length vector to implement spherical principle surface
    focConst = F**2+Rin**2
    focVec = sqrt(focConst-R**2)

    #Weight vector for shell radii
    weights = np.zeros(M*N)
    spanv = np.zeros(M)
    for i in range(M):
        #Full angular span of each shell
        spanv[i] = 2*np.arcsin(span/2/R[i])
        #Geometric area in each shell
        weights[i*N:(i+1)*N] = ((R[i]+d)**2-R[i]**2) * spanv[i]/2 / 10000.

    #Trace rays through SPO modules
    rays,refl = traceSPO(R,L,focVec,N,M,spanv,wave,d=d,t=t)
    weights = weights*refl

    #Determine outermost radius of grating array
    PT.transform(rays,0,0,focVec[-1]-(L.max()+gap+95.),0,0,0)
    PT.flat(rays)
    outerrad = np.max(sqrt(rays[1]**2+rays[2]**2))
    hubdist = sqrt(outerrad**2 + (focVec[-1]-(L.max()+gap+95.))**2)
    angle = np.arctan(outerrad/(focVec[-1]-(L.max()+gap+95.)))
    thetag = angle - 1.5*pi/180.
    print 'Outerrad: %f\nHubdist: %f\nLmax: %f\nOuter Focus: %f\n' % \
          (outerrad,hubdist,L.max(),focVec[-1])

    #Trace grating array - need to add in grating efficiency and
    #CCD quantum efficiency
    if bestFocus is None:
        return gratArray(rays,outerrad,hubdist,angle,inc,l=l,\
                         weights=weights,order=-3,blazeYaw=blazeYaw,\
                         wave=2.4)
##        return traceSector(Rin,Rout,F,N,span=span,d=d,t=t,gap=gap,\
##                inc=inc,l=l,bestFocus=bestFocus,order=order,\
##                blazeYaw=blazeYaw,wave=wave,marg=marg)
        
    gratArray(rays,outerrad,hubdist,angle,inc,l=l,bestFocus=bestFocus,\
                  weights=weights,order=order,blazeYaw=blazeYaw,wave=wave)

    #Account for grating efficiency
    geff = gratEff(order)
    #Add in grating efficiency and CCD QE
    weights = weights * geff(wave) * ccdQE(wave)

    #Go to focal plane
    tran.transform(rays,0,0,bestFocus,0,0,0)
    surf.flat(rays)

    #Get rid of rays that made it through
##    ind = rays[1] > 0
##    rays = PT.vignette(rays,ind=ind)
##    weights = weights[ind]
    #Get rid of outliers
    ind = np.abs(rays[2]-np.average(rays[2]))<10.
    rays = PT.vignette(rays,ind=ind)
    weights = weights[ind]

    #If no rays made it through
    if np.size(rays[1])==0:
        return 0,0,0

    #Look at resolution of this line
    try:
        cy = np.average(rays[2],weights=weights)
    except:
        pdb.set_trace()
    if findMargin is True:
        #Compute FWHM after Gaussian convolution
        fwhms = np.linspace(1,3,201)/60.**2*pi/180*12e3
        tot = np.array([convolveLSF(rays,.001,m,weights=weights)\
                        for m in fwhms/2.35])
        #Convert to arcsec, return std of margin convolution
        marg = fwhms[np.argmin(np.abs(tot/12e3*180/pi*60**2-3.))]/2.35
        return marg,rays
    fwhm = convolveLSF(rays,.001,marg,weights=weights)
    resolution = cy/fwhm
    area = np.sum(weights)*.799*.83*.8*.8#Grat plates, azimuthal ribs, packing^2
    #print resolution
    #print area
    #print sqrt((cy/3000)**2 - lsf**2)/F * 180/pi*60**2
    
    return resolution, area, np.nanmean(rays[1]), np.nanmean(rays[2]),\
           rays,weights

def traceSPO(R,L,focVec,N,M,spanv,wave,d=.605,t=.775):
    """Trace SPO surfaces sequentially. Collect rays from
    each SPO shell and set them to the PT rays at the end.
    Start at the inner radius, use the wafer and pore thicknesses
    to vignette and compute the next radius, loop while
    radius is less than Rout.
    """
    #Ray bookkeeping arrays
    trays = [np.zeros(M*N) for n in range(10)]

    #Loop through shell radii and collect rays
    ref = np.zeros(M*N)
    for i in range(M):
        #Set up source annulus
        rays = PT.subannulus(R[i],R[i]+d,spanv[i],N)
        z,n = rays[3],rays[6]
        #Transform rays to be above xy plane
        PT.transform(rays,0,0,0,pi,0,0)#n = -n
        PT.transform(rays,0,0,-100.,0,0,0)
        #Trace to primary
        PT.spoPrimary(rays,R[i],focVec[i])
        PT.reflect(rays)
        #Compute reflectivity
        inc = PT.grazeAngle(rays)#np.arcsin(l*ux+m*uy+n*uz)
        refl = CXCreflIr(inc,1239.8/wave,.5)
        #Vignette
        ind  = np.logical_and(rays[3]<=L[i],rays[3]>=0.)
        if np.sum(ind) < N:
            pdb.set_trace()
        PT.vignette(rays,ind)
        #Trace to secondary
        PT.spoSecondary(rays,R[i],focVec[i])
        PT.reflect(rays)
        #Compute reflectivity
        inc = PT.grazeAngle(rays)#inc = np.arcsin(l*ux+m*uy+n*uz)
        ref[i*N:(i+1)*N] = refl * CXCreflIr(inc,1239.8/wave,.5)
        #Vignette
        ind  = np.logical_and(rays[3]<=0.,rays[3]>=-L[i])
        if np.sum(ind) < N:
            pdb.set_trace()
        PT.vignette(rays,ind)
        #Set plane to be at focus
        PT.transform(rays,0,0,-focVec[i],0,0,0)
        #Collect rays
        try:
##            tx[i*N:(i+1)*N] = PT.x
##            ty[i*N:(i+1)*N] = PT.y
##            tz[i*N:(i+1)*N] = PT.z
##            tl[i*N:(i+1)*N] = PT.l
##            tm[i*N:(i+1)*N] = PT.m
##            tn[i*N:(i+1)*N] = PT.n
##            tux[i*N:(i+1)*N] = PT.ux
##            tuy[i*N:(i+1)*N] = PT.uy
##            tuz[i*N:(i+1)*N] = PT.uz
            for t in range(1,7):
                temp = trays[t]
                temp[i*N:(i+1)*N] = rays[t]
        except:
            pdb.set_trace()

    return trays,ref

def gratArray(rays,outerrad,hubdist,angle,inc,l=95.,bestFocus=None,\
              weights=None,order=0,blazeYaw=0.,wave=1.):
    """Trace rays leaving SPO petal to the fanned grating array.
    Start with outermost radius and rotate grating array about
    the hub. Define outermost grating position by max ray radius
    at desired axial height.
    Rays have been traced to bottom of outermost grating.
    """
    x,y = rays[1:3]
    #Dummy rays to ensure return of reference frame
##    rays2 = sources.subannulus(220.,223.,10.*pi/180,100)
    #Put origin at bottom of outermost grating
    PT.transform(rays,outerrad,0,0,0,0,0)
##    PT.transform(rays2,outerrad,0,0,0,0,0)
    #Go to proper incidence angle of grating
    PT.transform(rays,0,0,0,0,0,-pi/2)
    PT.transform(rays,0,0,0,-pi/2-angle+inc,0,0)
##    PT.transform(rays2,0,0,0,0,0,-pi/2)
##    PT.transform(rays2,0,0,0,-pi/2-angle+inc,0,0)
    #Go to hub
    PT.transform(rays,0,0,0,0,0,blazeYaw) #Put in blaze
    PT.transform(rays,0,hubdist,0,0,0,0)
##    PT.transform(rays2,0,0,0,0,0,blazeYaw) #Put in blaze
##    PT.transform(rays2,0,hubdist,0,0,0,0)
    #Trace out gratings until no rays hit a grating
    #Flat
    #Indices
    #Reflect
    #Apply Grating
    #Next
    PT.flat(rays)
    rho = -sqrt(x**2+y**2)*np.sign(y)
    ind = np.logical_and(rho>hubdist,rho<l+hubdist)
    ind2 = np.copy(ind)
    ang = l*sin(inc)/hubdist*.95
    
    i = 0
    prev = np.copy(ind)
    while np.sum(ind2)>0:
        i = i+1
##        pdb.set_trace()
        PT.reflect(rays,ind=ind2)            
        tran.radgrat(rays,160./hubdist,order,wave,ind=ind2)
        PT.transform(rays,0,0,0,ang,0,0)
##        PT.transform(rays2,0,0,0,ang,0,0)
        PT.flat(rays)
        rho = -sqrt(x**2+y**2)*np.sign(y)
        prev = np.logical_or(prev,ind) #Add rays hitting new grating
        ind = np.logical_and(rho>hubdist,rho<l+hubdist)
        #ind = np.logical_and(PT.y<-hubdist,PT.y>-l-hubdist)
        ind2 = np.logical_and(np.invert(prev),ind) #Remove previous rays
        #sys.stdout.write('%i \r' % i)
        #sys.stdout.flush()

##    #Go to focal plane
##    PT.transform(rays,0,-hubdist,0,0,0,0)
##    PT.transform(rays,0,0,0,0,0,-blazeYaw) #Reverse blaze
##    #Currently at bottom point of innermost grating

    #Get back to original outermost grating reference frame
    PT.transform(rays,0,0,0,-ang*i,0,0)
    PT.transform(rays,0,-hubdist,0,0,0,0)
    PT.transform(rays,0,0,0,0,0,-blazeYaw)
    PT.transform(rays,0,0,0,pi/2+angle-inc,0,0)
    PT.transform(rays,0,0,0,0,0,pi/2)
    PT.transform(rays,-outerrad,0,0,0,0,0)

##    PT.transform(rays2,0,0,0,-ang*i,0,0)
##    PT.transform(rays2,0,-hubdist,0,0,0,0)
##    PT.transform(rays2,0,0,0,0,0,-blazeYaw)
##    PT.transform(rays2,0,0,0,pi/2+angle-inc,0,0)
##    PT.transform(rays2,0,0,0,0,0,pi/2)
##    PT.transform(rays2,-outerrad,0,0,0,0,0)
    #Should be there
    PT.flat(rays)
    
##    PT.transform(rays,0,hubdist,0,0,0,0)
##    PT.transform(rays,0,0,0,-ang*i+pi/2+angle-inc,0,0)
##    PT.transform(rays,0,0,0,0,0,pi/2)
##    PT.flat(rays)
    
    #Find focus
    if bestFocus is None:
        return surf.focusY(rays,weights=weights)
##
##    #Focus already found, tracing diffracted line
##    PT.transform(rays,0,0,bestFocus,0,0,0)
##    PT.flat(rays)

    return None
    

def arc(inc,yaw,hubdist,wave,order,dpermm):
    """Return x and y positions of diffraction arc
    as a function of wavelength for a given order"""
    #Set up source ray
    PT.circularbeam(0.,1)
    #Transform to grating frame
    PT.transform(0,0,0,pi/2+inc,0,0)
    PT.transform(0,0,0,0,0,yaw)
    PT.flat()
    #Apply grating
    PT.reflect()
    PT.radgrat(hubdist,dpermm,order,wave)
    #Go to focus
    PT.transform(0,0,0,0,0,-yaw)
    PT.transform(0,hubdist,0,0,0,0)
    PT.transform(0,0,0,pi/2,0,0)
    PT.flat()
    #Get ray positions
    return PT.x[0],PT.y[0]

def plotArc(inc,yaw,hubdist,dpermm,order):
    """Plot arc vs. x and wavelength in a double plot"""
    #Define wavelength vector
    wave = np.linspace(.1,10.,1000)
    #Figure out coordinates for each order
    x = np.zeros(np.size(wave))
    y = np.copy(x)
    for i in range(np.size(wave)):
        x[i],y[i] = arc(inc,yaw,hubdist,wave[i],-1,dpermm)

    pdb.set_trace()
    #Wavelength scale
    scale = np.nanmean(wave/x/order)
    fn = lambda x: scale*x
    plotting.pltd2(x,y,fn)

    return x,y,scale

def effPlot(wave,w1,w2):
    """Make grating efficiency plot
    Indicate order ranges using low 1st order wave w1
    and high 1st order wave w2
    """
    color = ['b','r','g','c','k','y','m']
    for order in np.arange(-7,-1):
        ef = gratEff(order)
        plt.plot(wave,ef(wave),label=str(order),color=color[order+7])
        #plt.plot([-w1/order,-w1/order],[0.,.7],'--',color=color[order+7])
        #plt.plot([-w2/order,-w2/order],[0.,.7],'--',color=color[order+7])

def predictionPlots(wave,arcus,w1,w2):
    """Make plots of resolution and effective area
    Plot only regions covered by CCDs
    """
    eafig = plt.figure()
    plt.grid()
    resfig = plt.figure()
    plt.grid()
    for i in range(np.shape(arcus)[0]):
        order = i+1
        ind = np.logical_and(wave<=w2/order,wave>=w1/order)
        plt.figure(resfig.number)
        plt.plot(wave[ind],arcus[i,0][ind],label=str(order))
        plt.figure(eafig.number)
        plt.plot(wave[ind],arcus[i,1][ind]*4,label=str(order))

def plotPetal(x,y,tra=tr.identity_matrix(),color='blue'):
    """Make an isometric plot of a petal layout"""
    #Apply transformation matrices
    arc = [x,y,np.zeros(len(x)),np.ones(len(x))]
    arc = np.dot(tra,arc)
    #Plot arc
    plotting.isoplot(arc[0],arc[1],color=color)
    #Plot petal aperture
    appx = np.array([300.,300.,800.,800.,300.])
    appy = np.array([187.5,-187.5,-187.5,187.5,187.5])
    app = [appx,appy,np.zeros(5),np.ones(5)]
    app = np.dot(tra,app)
    plt.plot(app[0],app[1],color=color)
    #Plot zero order
    zo = np.dot(tra,[615.253,0,0,1])
    plt.plot(zo[0],zo[1],'.',color=color)
    return [arc,app,zo]

def matchArc(arc):
    """Take the raytraced arc, negate the dispersion direction, and
    compute the translation and rotation transformations to
    match the two spectra
    """
    #Create the second arc by negating the dispersion direction
    sh = np.shape(arc)
    arc = [arc[0],arc[1],np.zeros(sh[1]),np.ones(sh[1])]
    arc2 = np.copy(arc)
    arc2[1] = -arc2[1]

    #Transform second arc by pi/2
    r = tr.rotation_matrix(pi/2,[0,0,1])
    arc2 = np.dot(r,arc2)

    #Get gradients
    dx = np.diff(arc2[0])[0]
    gr = np.gradient(arc2[1],dx)

    #Apply rotation to match slopes
    angle = pi/2+np.arctan(gr[-1])+np.arctan(gr[0])
    r2 = tr.rotation_matrix(-angle,[0,0,1])
    arc3 = np.dot(r2,arc2)

    #Apply translation to match 
    pdb.set_trace()
##    arc3[0] = arc3[0] + (arc[0][-1]-arc3[0][0])
##    arc3[1] = arc3[1] + (arc[1][-1]-arc3[1][0])
    t = tr.translation_matrix([(arc[0][-1]-arc3[0][0]),\
                               (arc[1][-1]-arc3[1][0]),\
                               0.])
    arc3 = np.dot(t,arc3)
    tra = np.dot(t,(np.dot(r2,r)))

    #Get even closer
    dist = np.sqrt((arc[0][0]-arc3[0][-1])**2+(arc[1][0]-arc3[1][-1])**2)
    lever = np.sqrt((arc[0][0]-arc[0][-1])**2+(arc[1][0]-arc[1][-1])**2)
    angle = dist/lever
    tra2 = np.dot(tr.rotation_matrix(angle,[0,0,1],point=\
                                     [arc3[0][0],arc3[1][0],0]),tra)

    #Return transformation
    return tra2,angle

def makeLayout(arc,sepAngle):
    tra,angle = matchArc(arc) #Get nominal transformation for paired petal
    tra2 = np.dot(tr.rotation_matrix(angle/2,[0,0,1]),tra)
    tra1 = tr.rotation_matrix(angle/2,[0,0,1])

    #Create the second arc by negating the dispersion direction
    sh = np.shape(arc)
    arc = [arc[0],arc[1],np.zeros(sh[1]),np.ones(sh[1])]
    arc2 = np.copy(arc)
    arc2[1] = -arc2[1]

    #Apply transformations
    arc2 = np.dot(tra2,arc2)
    arc1 = np.dot(tra1,arc)
    
    #Find center of second arc
    xc,yc = fit.circle(arc2[0],arc2[1],0,0)[0]

    #Rotate each petal about this point
    tra1 = np.dot(tr.rotation_matrix(sepAngle,[0,0,1],point=[xc,yc,0]),tra1)
    tra2 = np.dot(tr.rotation_matrix(-sepAngle,[0,0,1],point=[xc,yc,0]),tra2)

    #Add translation to second petal
    tra2 = np.dot(tr.translation_matrix([10/np.sqrt(2),10/np.sqrt(2),0]),tra2)

    #Now plot the petals with the correct translations
    pet1 = plotPetal(arc[0],arc[1],tra=tra1,color='blue')
    pet2 = plotPetal(arc[0],-arc[1],tra=tra2,color='red')

    #Plot dashed line between corners closest to arcs
    x1,x2,y1,y2=pet1[1][0][3],pet2[1][0][2],pet1[1][1][3],pet2[1][1][2]
    plt.plot([pet1[1][0][3],pet2[1][0][2]],[pet1[1][1][3],pet2[1][1][2]],'k--')

def convolveLSF(rays,binsize,std,weights=None,plot=False):
    """Convolve a gaussian with the LSF determined by
    histogramming the ray LSF
    std should be supplied as a distance in mm"""
    #Bin up rays in dispersion direction
    n,b = np.histogram(rays[2],bins=\
                       np.arange(rays[2].mean()-.5,rays[2].mean()+.5,binsize),\
                       weights=weights)
    b = np.array([np.mean([b[i],b[i+1]]) for i in range(len(b)-1)])
    #Create convolution kernel
    gaussk = conv.Gaussian1DKernel(std/binsize)
    n2 = conv.convolve(n,gaussk)
    #Determine FWHM
    maxi = np.argmax(n2) #Index of maximum value
    #Find positive bound
    bp = b[maxi:]
    fwhmp = bp[np.argmin(np.abs(n2[maxi:]-n2.max()/2))]-b[maxi]
    bm = b[:maxi]
    fwhmm = b[maxi]-bm[np.argmin(np.abs(n2[:maxi]-n2.max()/2))]
    if plot is True:
        plt.plot(b,n)
        plt.plot(b,n2)
    return fwhmm+fwhmp
