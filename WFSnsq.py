#Script to analyze detector data from Zemax non-sequential optics
#WFS simulation.  Need to pull out list of centroid locations to reduce
#data size.  This is similar to CCD 'photon counting' algorithm from
#the GSFC resolution test

from numpy import *
from matplotlib.pyplot import *
import pdb, sys, pickle, struct
import reconstruct as frec

def centroidlist(filename,thresh):
    #Read in Zemax detector data
    d = transpose(genfromtxt(filename,skip_header=24))
    d = d[1:]

    #Find maximum pixel, calculate centroid of pixels surrounding
    #that pixel, save centroid coordinates and set pixels to zero
    xc = []
    yc = []
    while d.max()>thresh:
        #Find maximum
        xi,yi = where(d==d.max()) #Coordinates of max pixels
        xi = xi[0]
        yi = yi[0] #Handles multiple maxima - take first

        #Calculate centroid in box +- 3 pixels from max
        #Create X and Y arrays to simplify calculation
        x,y = meshgrid(arange(xi-3,xi+4),arange(yi-3,yi+4))
        dsub = d[xi-3:xi+4,yi-3:yi+4]
        xc.append(sum(x*dsub)/sum(dsub))
        yc.append(sum(y*dsub)/sum(dsub))

        #Set pixels to zero
        d[xi-3:xi+4,yi-3:yi+4] = 0.

        sys.stdout.write(str(size(xc)/128.**2) + '\n')

    #Return centroid list
    return xc, yc

#Function to read in Zemax ray database file in text format
def readrays(filename):
    #Open file in read mode
    f = open(filename,'r')

    #Loop through file, determine if line is of interest, and
    #save position and direction cosines
    x = []
    y = []
    l = []
    m = []
    n = []
    li = f.readline()
    while li != '':
        #Is line full segment (split size greater than 17)
        lis = li.split()
        if size(lis) > 16:
            #Is line intercept at detector (hits object 7)
            if lis[4] == '8':
                #Save x,y coordinates
                x.append(lis[9])
                y.append(lis[10])
                #Save l,m,n cosines
                l.append(lis[12])
                m.append(lis[13])
                n.append(lis[14])
        #Read in next line
        li = f.readline()
        sys.stdout.write(li+'\r')

    #Return list of coordinates and direction cosines
    return array(x).astype('double'),array(y).astype('double'),\
           array(l).astype('double'),array(m).astype('double'),\
           array(n).astype('double')

#Bin up ray intercepts into lenslet pixels and output array of average
#direction cosines in x and in y for each lenslet
def WFSbin(filename):
    #Read in rays
    x,y,l,m,n = readrays(filename)

    #What is radius of image?
    rad = sqrt(max(x**2+y**2))

    #Set up lenslet binning, simulate 128 by 128 HASO
    #Want 2d arrays for xslope and yslope
    #Also save number of rays per lenslet
    xslope = zeros((128,128))
    yslope = copy(xslope)
    xstd = copy(xslope)
    ystd = copy(yslope)
    n = copy(xslope)
    for xi in range(128):
        xpos = (xi - 64) * .114 + .114/2 #X position of lenslet
        xsel = logical_and(x>xpos-.114/2,x<xpos+.114/2) #X slice
        for yi in range(128):
            ypos = (yi - 64) * .114 + .114/2 #Y position of lenslet
            ysel = logical_and(y>ypos-.114/2,y<ypos+.114/2) #Y slice
            #Select rays in lenslet and calculate xslope,yslope
            ind = logical_and(xsel,ysel)
            if sum(ind) > 0:
                xslope[xi,yi] = mean(l[ind])
                yslope[xi,yi] = mean(m[ind])
                n[xi,yi] = sum(ind)
                xstd[xi,yi] = std(l[ind])
                ystd[xi,yi] = std(m[ind])

    #Return lenslet data
    return xslope, yslope, n, xstd, ystd

#Define class for WFS results
class wave:
    def __init__(self,xs,ys,n):
        self.xs = xs
        self.ys = ys
        self.n = n

#Load in waves and analyze
def analyzewaves():
    #Load raytraced waves
    f = open('140210FlatWave.txt','r')
    flat = pickle.load(f)
    f.close()
    f = open('140210FlatInf.txt','r')
    flatinf = pickle.load(f)
    f.close()
    f = open('140210DeformedWave.txt','r')
    deformed = pickle.load(f)
    f.close()
    f = open('140210DeformedInf.txt','r')
    deformedinf = pickle.load(f)
    f.close()

    #Compute differences in slopes at each lenslet
    #Only take lenslets with non-zero slopes
    #Set zero slopes to NaNs
    flatdiff = zeros((128,128))
    deformeddiff = zeros((128,128))
    for xi in range(128):
        for yi in range(128):
            if flat.xs[xi,yi]==0:
                flat.xs[xi,yi]=NaN
                flat.ys[xi,yi]=NaN
            if flatinf.xs[xi,yi]==0:
                flatinf.xs[xi,yi]=NaN
                flatinf.ys[xi,yi]=NaN
            if deformed.xs[xi,yi]==0:
                deformed.xs[xi,yi]=NaN
                deformed.ys[xi,yi]=NaN
            if deformedinf.xs[xi,yi]==0:
                deformedinf.xs[xi,yi]=NaN
                deformedinf.ys[xi,yi]=NaN

    #Save NaN arrays to disk
    f = open('140210FlatWave.txt','w')
    pickle.dump(flat,f)
    f.close()
    f = open('140210FlatInf.txt','w')
    pickle.dump(flatinf,f)
    f.close()
    f = open('140210DeformedWave.txt','w')
    pickle.dump(deformed,f)
    f.close()
    f = open('140210DeformedInf.txt','w')
    pickle.dump(deformedinf,f)
    f.close()

    flatdiffx = flat.xs/flat.n-flatinf.xs/flatinf.n
    flatdiffy = flat.ys/flat.n-flatinf.ys/flatinf.n
    defdiffx = deformed.xs/deformed.n-deformedinf.xs/deformedinf.n
    defdiffy = deformed.ys/deformed.n-deformedinf.ys/deformedinf.n

    #Plot differences
    figure(1)
    clf()
    contourf(flatdiffx,levels=linspace(-.0001,.0001\
                                       ,100))
    figure(2)
    clf()
    contourf(defdiffx,levels=linspace(-.0001,.0001\
                                      ,100))

    figure(3)
    clf()
    diffdiffx = flatdiffx-defdiffx
    diffdiffy = flatdiffy-defdiffy
    contourf(diffdiffx,levels=linspace(-.0001,.0001,100))
    figure(4)
    clf()
    contourf(diffdiffy,levels=linspace(-.0001,.0001,100))
    
    
#Read in ZRD binary file
def readZRD(filename,xc=0.,yc=0.):
    f = open(filename,'rb')

    #Skip header data
    f.seek(8)

    #Read in number of segments in ray, if ='', end of file
    ray = 1
    xang = zeros((128,128))
    yang = copy(xang)
    n = copy(xang)
    segbytes = f.read(4)
    while segbytes != '':
        numseg = struct.unpack('i',segbytes)[0] #Number of segments
        for i in range(numseg):
            f.seek(8,1) #Skip to hit object
            obj = struct.unpack('i',f.read(4))[0] #Which object was hit?
            if obj != 8:
                f.seek(7*4+21*8,1) #Skip to next segment
            else:
                #Read in ray data
                f.seek(7*4+2*8,1) #Skip to x position
                x = struct.unpack('d',f.read(8))[0]-xc #Read x
                y = struct.unpack('d',f.read(8))[0]-yc #Read y
                f.seek(8,1) #Skip to cosines
                l = struct.unpack('d',f.read(8))[0] #Read x cosine
                m = struct.unpack('d',f.read(8))[0] #Read y cosine
                #Put ray into proper bin
                xbin = floor(x/.114) + 65 #X bin
                ybin = floor(y/.114) + 65 #Y bin
                
                #Add to cumulative arrays
                n[xbin,ybin] += 1
                xang[xbin,ybin] += l
                yang[xbin,ybin] += m

                f.seek(14*8 + 208*(numseg-i-1),1) #Skip to next ray
                segbytes = f.read(4) #Read number of segments in next ray

                sys.stdout.write('Ray ' + str(ray) + '\r') #Print status
                sys.stdout.flush()
                ray += 1
                break #Break out of loop

    nanmask = where(n==0)
    xang[nanmask] = nan
    yang[nanmask] = nan
    xang = tan(arcsin(xang / n))
    yang = tan(arcsin(yang / n))

    return wave(xang,yang,n)

#Phase reconstruction algorithm (Southwell 1980)
def reconstruct(wv,criteria):
    #Initialize phase array to all NaNs
    phase = zeros((130,130))

    #Initialize slope arrays
    xang = reshape(repeat(NaN,130*130),(130,130))
    yang = copy(xang)
    for xi in arange(1,129):
        for yi in arange(1,129):
            xang[xi,yi] = wv.xs[xi-1,yi-1]
            yang[xi,yi] = wv.ys[xi-1,yi-1]
    nanmask = isnan(xang)

    #Set non-existent phase elements to NaN
    phase[nanmask] = nan
    phasec = copy(phase)

    h = 114. #114 um step size

    w = 2/(1+sin(pi/(128+1))) #Relaxation factor

    #Enter reconstruction loop, only loop through non-nan elements
    rms = 1.
    while rms > criteria: #Convergence condition is 100th wave RMS improvement
        #Loop through phase array and update in Gauss-Seidel fashion
        for xi in arange(1,129):
            for yi in arange(1,129):
                #If NaN, continue to next element
                if isnan(phasec[xi,yi]):
                    continue
                
                #Make array of nearest neighbors
                xneigh = array([xang[xi+1,yi],xang[xi-1,yi],\
                                xang[xi,yi+1],xang[xi,yi-1]])
                yneigh = array([yang[xi+1,yi],yang[xi-1,yi],\
                                yang[xi,yi+1],yang[xi,yi-1]])
                nanneighb = isnan(xneigh)
                xneigh[nanneighb] = -xang[xi,yi]
                yneigh[nanneighb] = -yang[xi,yi]
                
                #Compute slope term
                bk = .5*(yneigh[2]-yneigh[3]+\
                         xneigh[0]-xneigh[1])*h
                
                #Compute nearest neighbor sum
                psum = nansum(array([phasec[xi+1,yi],phasec[xi-1,yi],\
                           phasec[xi,yi+1],phasec[xi,yi-1]]))
                
                #If lenslet has no good neighbors, set to NaN to discard
                if isnan(psum):
                    phasec[xi,yi] = nan
                    phase[xi,yi] = nan
                    xang[xi,yi] = nan
                    yang[xi,yi] = nan
                    continue
                
                #Compute gk factor
                if (xi==1 or xi==128) and (yi==1 or yi==128):
                    gk = 2
                elif (xi==1 or xi==128) and (1<yi<128):
                    gk = 3
                elif (yi==1 or yi==128) and (1<xi<128):
                    gk = 3
                else:
                    gk = 4
                #Correct phase
                phasec[xi,yi] = phasec[xi,yi] + w*(psum/gk+bk/gk-phasec[xi,yi])
                
        #Determine rms improvement
        notnan = invert(isnan(xang))
        rms = sqrt(mean((phasec[notnan]-phase[notnan])**2))
##        rms = abs(phasec[notnan]-phase[notnan])
##        rms = rms.max()
        sys.stdout.write('RMS: ' + str(rms))
        sys.stdout.flush()
        phase = copy(phasec)

    return phasec

#Phase reconstruction algorithm (Southwell 1980)
#Calls reconstruct algorithm written in Fortran
def freconstruct(wv,criteria):
    #Initialize phase array to all NaNs
    phase = reshape(repeat(nan,130*130),(130,130))

    #Initialize slope arrays
    xang = reshape(repeat(100.,130*130),(130,130))
    yang = copy(xang)
    for xi in arange(1,129):
        for yi in arange(1,129):
            xang[xi,yi] = wv.xs[xi-1,yi-1]
            yang[xi,yi] = wv.ys[xi-1,yi-1]
            phase[xi,yi] = 0.
    nanmask = isnan(xang)

    #Reset nans to 100
    xang[nanmask] = 100.
    yang[nanmask] = 100.

    #Set non-existent phase elements to NaN
    phase[where(xang==100)] = 100.

    pdb.set_trace()

    #Enter reconstruction loop, only loop through non-nan elements
    phasec = frec.reconstruct(xang,yang,criteria,phase)

    #Turn padded 0's into nans
    phasec[where(phasec==100.)] = nan

    return phasec

#Create wavefront arrays from ray databases
def parsetraces():
    xang,yang,n = readZRD('140213IdealFlatTrace.zrd')
    xang[where(n==0)] == nan
    yang[where(n==0)] == nan
    xang = xang/n
    yang = yang/n
    f = open('140213IdealFlatSlopes.txt','w')
    pickle.dump(wave(xang,yang,n),f)
    f.close()
    xang,yang,n = readZRD('140213IdealInfluenceTrace.zrd')
    xang[where(n==0)] == nan
    yang[where(n==0)] == nan
    xang = xang/n
    yang = yang/n
    f = open('140213IdealInfluenceSlopes.txt','w')
    pickle.dump(wave(xang,yang,n),f)
    f.close()
    xang,yang,n = readZRD('140213RNR283Trace.zrd')
    xang[where(n==0)] == nan
    yang[where(n==0)] == nan    
    xang = xang/n
    yang = yang/n
    f = open('140213RNR283Slopes.txt','w')
    pickle.dump(wave(xang,yang,n),f)
    f.close()
    xang,yang,n = readZRD('140213RNR283InfluenceTrace.zrd')
    xang[where(n==0)] == nan
    yang[where(n==0)] == nan
    xang = xang/n
    yang = yang/n
    f = open('140213RNR283InfluenceSlopes.txt','w')
    pickle.dump(wave(xang,yang,n),f)
    f.close()

#Reconstruct wavefront in reference mode
def refconstruct(meas,ref):
    inf = wave(meas.xs-ref.xs,meas.ys-ref.ys,meas.n+ref.n)
    return freconstruct(inf,1.e-12)

#Construct difference slope files and reconstruct influence functions
def influence():
    #Create ideal flat slopes
    f = open('140213IdealFlatSlopes.txt','r')
    flat = pickle.load(f)
    f.close()
    f = open('140213IdealInfluenceSlopes.txt','r')
    flatinf = pickle.load(f)
    f.close()
    f = open('140213RNR283Slopes.txt','r')
    deformed = pickle.load(f)
    f.close()
    f = open('140213RNR283InfluenceSlopes.txt','r')
    deformedinf = pickle.load(f)
    f.close()

    flatwave = freconstruct(flat,1.e-12)
    flatact = freconstruct(flatinf,1.e-12)
    rnrwave = freconstruct(deformed,1.e-12)
    rnract = freconstruct(deformedinf,1.e-12)

    flatslopes = wave(flatinf.xs-flat.xs,flatinf.ys-flat.ys,flatinf.n+flat.n)
    defslopes = wave(deformedinf.xs-deformed.xs,deformedinf.ys-deformed.ys,\
                     deformedinf.n+deformed.n)
    rnr283slopes = wave(deformed.xs-flat.xs,deformed.ys-flat.ys,flat.n+deformed.n)

    #Reconstruct wavefronts
    flatinf = freconstruct(flatslopes,1.e-12)
    definf = freconstruct(defslopes,1.e-12)
    rnr283 = freconstruct(rnr283slopes,1.e-12)

    return flatinf, definf, flatwave, flatact, rnrwave, rnract, rnr283

#Code to determine image location and magnification of 20 mm beam system
def sys20mm(so1=100.,d=400.,f1=400.,f2=200.):
    #Input lens
    s01 = 100. #Object 10 cm from collimator

    #Compute first image location and magnification
    si1 = (1/f1-1/so1)**(-1)
    mag1 = -si1/so1
    print si1,mag1

    #Feed image into eyepiece
    so2 = d - si1
    si2 = (1/f2-1/so2)**(-1)

    #Total magnification
    mag2 = -si2/so2
    magt = mag2*mag1
    print si2/so2

    return si2,magt
    
