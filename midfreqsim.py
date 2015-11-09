import traces.PyTrace as PT
from numpy import *
from matplotlib.pyplot import *
import pdb,time,reconstruct
from utilities.plotting import mycontour,nanmean
import zernikemod as zmod
import traces.conicsolve as conicsolve

#Trace rays from test optic plane through lens system to WFS
def tracefromtest(fieldx=0.,fieldy=0.,fieldz=0.,imgx=0.,imgy=0.,imgz=0.,\
                  colx=0.,coly=0.,colz=0.):
    #Trace to collimator lens
    PT.transform(0.,0.,100.,0.,0.,0.)
    PT.transform(colx,coly,colz,0,0,0)
    PT.flat() #Get to center of first surface
    PT.lens(206.03,0.,4.57,25.4*2,1.51501) #First lens
    PT.transform(-colx,-coly,-colz,0,0,0)

    #Trace to beamsplitter
    PT.transform(0.,0.,283.602,0,0,0)
    PT.flat()
    PT.refract(1.,1.51501)#Refract into beamsplitter
    PT.transform(0,0,50,0,0,0)
    PT.flat()
    PT.refract(1.51501,1.)#Refract out of beamsplitter

    #Trace to field lens
    PT.transform(0,0,211.828,0,0,0)
    PT.transform(fieldx,fieldy,fieldz,0,0,0) #Field lens misalignment
    PT.flat()
    PT.lens(205.86,-205.86,5.05,50.,1.51501)
    PT.transform(-fieldx,-fieldy,-fieldz,0,0,0) #Reverse lens misalignment

    #Trace to image plane
    PT.transform(0.,0.,285.23482227,0,0,0)
    PT.transform(imgx,imgy,imgz,0,0,0) #Include WFS misalignment
    PT.flat()

#Set up central point source for imaging check
def centralpointsource(num,testx=0.,testy=0.,testz=0.,**kwgs):
    PT.pointsource(pi/180,num)
    PT.transform(-testx,-testy,-testz,0,0,0)
    tracefromtest(**kwgs)

#Define circular beam and bounce off reference flat
def reference(num,pitch=0.,yaw=0.,roll=0.,**kwgs):
    PT.circularbeam(12.5,num)
    PT.transform(0,0,0,pi,0,0)
    PT.transform(0,0,-100,0,0,0)

    PT.transform(0,0,0,pitch,yaw,roll)
    PT.flat()
    PT.reflect()
    PT.itransform(0,0,0,pitch,yaw,roll)

    tracefromtest(**kwgs)

#Test retrace error
def retrace(num,pitch=0.,yaw=0.,roll=0.,**kwgs):
    #Trace reference flat
    reference(10**6,**kwgs)
    #Move to center of image plane
    PT.x = PT.x - nanmean(PT.x)
    PT.y = PT.y - nanmean(PT.y)
    #Save reference slopes
    xang,yang,phase = reconstruct.southwellbin(PT.x,PT.y,PT.l,PT.m,.114,130,130)

    #Repeat with tilted test optic
    reference(10**6,pitch=pitch,yaw=yaw,roll=roll,**kwgs)
    #Move to center of image plane
    PT.x = PT.x - nanmean(PT.x)
    PT.y = PT.y - nanmean(PT.y)
    #Save reference slopes
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
    xanginf[ind] = 100
    yanginf[ind] = 100

    #Subtract average tip and tilt
    xanginf[invert(ind)] = xanginf[invert(ind)] - nanmean(xanginf[invert(ind)])
    yanginf[invert(ind)] = yanginf[invert(ind)] - nanmean(yanginf[invert(ind)])

    #Reconstruct influence wavefront
    influence = reconstruct.reconstruct(xanginf,yanginf,1.e-12,.114,phaseinf)

    #Make invalid pixels NaNs
    ind = where(influence==100)
    influence[ind] = NaN

    return influence

#Define annulus, trace perfect alignment, trace pitch and yaw
#Determine figure of merit based on edge shift
#RMS difference before and after misalignment
def edge(num,pitch=0.,yaw=0.,roll=0.,**kwgs):
    PT.edgebeam(12.5,10**4)
    PT.transform(0,0,0,pi,0,0)
    PT.transform(0,0,-100,0,0,0)

    PT.transform(0,0,0,pitch,yaw,roll)
    PT.flat()
    PT.reflect()
    PT.transform(0,0,0,-pitch,-yaw,-roll)
    
    tracefromtest(**kwgs)

def edgetest(num,pitch=0.,yaw=0.,roll=0.,**kwgs):
    #Reference
    edge(num,**kwgs)
    refx = PT.x
    refy = PT.y

    #Misalign pitch
    edge(num,pitch=pitch,yaw=yaw,**kwgs)

    #Figure of merit
    diff = sqrt(mean((refx-PT.x)**2+(refy-PT.y)**2))

    return diff

#Define circular beam and bounce off Wolter primary with sinusoid
def wolterripple(num,amp,freq,pitch=0.,yaw=0.,roll=0.,\
                 testx=0.,testy=0.,testz=0.,**kwgs):
    PT.rectbeam(2.,12.5,num)

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
    PT.woltersine(r0,z0,amp,freq)
    PT.reflect()
    #Go back to center of mirror and reverse misalignment
    PT.transform(0,-conicsolve.primrad(8476.,r0,z0),8476.\
                 ,-pitch,-yaw,-roll)
    #Rotate to nominal optical axis
    PT.transform(0,0,0,alpha-pi/2,0,0)
    #Put in spatial misalignment of test optic testz=-50 means
    #distance to collimator is 150 mm
    PT.transform(-testx,-testy,-testz,0,0,0)

    #Propagate through rest of system
    tracefromtest()

#Define circular beam and bounce off Wolter primary
def wolter(num,pitch=0.,yaw=0.,roll=0.,**kwgs):
    PT.circularbeam(12.5,num)

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

    #Propagate through rest of system
    tracefromtest()

#Reconstruct wavefront from rippled Wolter surface
def reconstructripple(amp,freq,testx=0.,testy=0.,testz=0.,**kwgs):
    #Trace reference flat
    reference(10**6,**kwgs)
    #Move to center of image plane
    PT.x = PT.x - nanmean(PT.x)
    PT.y = PT.y - nanmean(PT.y)
    #Save reference slopes
    xang,yang,phase = reconstruct.southwellbin(PT.x,PT.y,PT.l,PT.m,.114,130,130)

    #Trace rippled Wolter
    wolterripple(10**6,amp,freq,testx=testx,testy=testy,testz=testz,**kwgs)
    #Move to center of image plane
    PT.x = PT.x - nanmean(PT.x)
    PT.y = PT.y - nanmean(PT.y)
    #Get Wolter slopes
    xang2,yang2,phase2 = reconstruct.southwellbin(PT.x,PT.y,PT.l,PT.m,.114,130,130)

    #Construct phase and angle arrays for reconstruction of influence
    phaseinf = copy(phase)
    phaseinf[:,:] = 0.
    ind = where(logical_or(phase==100,phase2==100))
    phaseinf[ind] = 100
    xanginf = copy(xang)
    xanginf = xang2-xang
    yanginf = copy(yang)
    yanginf = yang2-yang
    xanginf[ind] = 100
    yanginf[ind] = 100

    #Reconstruct influence wavefront
    influence = reconstruct.reconstruct(xanginf,yanginf,1.e-12,phaseinf)

    #Make invalid pixels NaNs
    ind = where(influence==100)
    influence[ind] = NaN

    #Figure out pixel size
    centralslice = influence[65]
    centralslice = centralslice[invert(isnan(centralslice))]
    pixelsize = 25./size(centralslice)

    #Plot central psd - add loop to investigate other PSDs
    for i in range(shape(influence)[0]):
        l = influence[i]
        l = l[invert(isnan(l))]
        if size(l) > 50:
            f,s = axialPSD(l,pixelsize,window=False)
            if i % 3 == 0:
                plot(f,s,label=str(i))

    return influence,pixelsize

#Make PSD of an axial slice
#Dx just defines units of frequency
#Power is normalized to Parseval's theorem
#Divide by frequency interval when plotting
def axialPSD(l,dx,window=False):
    l = l[invert(isnan(l))]
    length = dx * (size(l)-1)
    x = arange(size(l))*dx
    fit = polyfit(x,l,2)
    l = l - polyval(fit,x)

    freq = fft.fftfreq(size(x),dx)
    if window is True:
        win = hanning(size(x))/sqrt(mean(hanning(size(x))**2))
    else:
        win = 1.
    spec = absolute(fft.fft(l*win))**2/size(l)**2
    spec = 2*spec[freq>0]
    freq = freq[freq>0]

    return freq,spec

#Determine ratio of image shift to optic roll
def comprollratio():
    wolter(10**4)
    refx = mean(PT.x)
    cx = []
    for roll in linspace(0.,2.,100):
        wolter(10**4,roll=roll*pi/180)
        cx.append(mean(PT.x)-refx)

    ratio = cx/linspace(0.,2.,100)
    return ratio

#Analyze image for yaw
def computeyaw(img,roll):
    #Loop through rows and find right edge
    edge = []
    for r in range(shape(img)[1]):
        #Make sure row is populated
        row = img[:,r]
        ind = where(row>0)[0]
        if size(ind)>5:
            if (roll < 0):
                #Find last non-zero pixel
                last = max(ind)
                norm = row[last-1] #This pixel is fully illuminated
                edge.append(last-1+row[last]/norm)
            if (roll > 0):
                #Find first non-zero pixel
                first = min(ind)
                norm = row[first+1]
                edge.append(first+1-row[first]/norm)

    x = arange(size(edge))
    fit = polyfit(x[3:-3],edge[3:-3],1)
    pdb.set_trace()

    return fit,edge

#Compute pitch, yaw, and roll of optic from WFS image
def computealign(num,pitch=0.,yaw=0.,roll=0.,**kwgs):
    #Trace reference flat
    reference(num,**kwgs)
    #Save image center
    cx = mean(PT.x)
    cy = mean(PT.y)
    #Save binned wavefront sensor data
    xang,yang,phase = reconstruct.southwellbin(PT.x,PT.y,PT.l,PT.m,.114,130,130)
    ind = where(phase==100.)
    xang[ind] = NaN
    yang[ind] = NaN
    phase[ind] = NaN

    #Trace Wolter primary with any misalignments
    wolter(num,pitch=pitch,yaw=yaw,roll=roll,**kwgs)
    #Compute roll by centroid in x direction
    xshift = mean(PT.x-cx)
    roll = xshift/2.07
    #Bin things up into an image
    img = hist2d(PT.x,PT.y,bins=arange(-7.3,7.3,.114))[0]
    #Analyze image to fit for yaw
    fit,edge = computeyaw(img,roll)
    yaw = -fit[0]*180/pi

    #Compute pitch by average y slope
    xang2,yang2,phase2 = reconstruct.southwellbin(PT.x,PT.y,PT.l,PT.m,.114,130,130)
    ind = where(phase2==100.)
    xang2[ind] = NaN
    yang2[ind] = NaN
    phase2[ind] = NaN
    #What index is central axial slice?
    cx,cy = zmod.locateimage(yang2,calcrad=False)
##    yslope = nanmean(yang2[round(cy)]-yang[round(cy)])
    yslope = nanmean(yang2-yang)
    #What is magnification? Compute max size in axial direction
    maxsize = 0
    for i in range(shape(img)[0]):
        axial = size(where(img[i]>0)[0])
        if (axial > maxsize):
            maxsize = axial
    mag = maxsize*.114/25.
    yslope = yslope*180/pi*mag/2 * 60
    pdb.set_trace()
    print yaw,roll,yslope
