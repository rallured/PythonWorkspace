from numpy import *
from matplotlib.pyplot import *
from plotting import mycontour,nanmean,pltd
import glob,os,pdb
from mpl_toolkits.mplot3d import Axes3D
import scipy.sparse as sparse
import sys
from zernikemod import *
from scipy.interpolate import griddata

#Need function to load in list of x,y,z
#positions for a given influence function
#First load in reference nodes
#Then offset them using the .dat files
def influence(num):
    direct = '/Users/ryanallured/GoogleDrive/WFS/SystemAlignment/RNR283/SquareFEA/'
    #Load in reference nodes
    nodes = transpose(genfromtxt(direct+'100mm_round_nodes.txt'))

    #Load in influence model
    f = glob.glob(direct+'*.dat')[num]
    inf = transpose(genfromtxt(f,skip_header=18,skip_footer=5))

    #Loop through nodes and apply offset
    for i in range(shape(nodes)[1]):
        for l in arange(1,4):
            nodes[l][i] = nodes[l][i] + inf[l][i]

    #Return x,y,z
    return nodes[1],nodes[2],nodes[3]

def influenceimage2(num,N):
    #Load in reference nodes
    x,y,z = influence(num)

    #Define inputs for griddata
    points = transpose((x,y))
    gridx,gridy = meshgrid(linspace(-50.,50.,N),linspace(-50.,50.,N))
    gridx = gridx.flatten()
    gridy = gridy.flatten()
    grid = griddata(points,z,(gridx,gridy),method='linear')

    return reshape(grid,(N,N))
    

#Load in influence function and convert to 2D image array
#Input influence function number
def influenceimage(num,N):
    xi,yi,zi = influence(num)

    #Create FEA model image array
    x,y = meshgrid(linspace(-50.,50.,N),\
                   linspace(-50.,50.,N))
    h = zeros(shape(x))
    for i in range(N):
        for j in range(N):
            #Find closest xi,yi position
            dist = (x[i,j]-xi)**2 + (y[i,j]-yi)**2
            if sqrt(x[i,j]**2+y[i,j]**2) < 50.:
                h[i,j] = zi[argmin(dist)]
            else:
                h[i,j] = NaN

    return h

#Construct summed influence function
def summedinf(N):
    inf = influenceimage2(0,N)
    for i in arange(48):
        inf = inf + influenceimage2(i+1,N)
        print i
        sys.stdout.flush()
    return inf

#Determine location of peak
def peaklocation(num):
    x,y,z = influence(num)

    figure()
    ax = axes(projection='3d')
    ax.plot3D(x,y,z)

    ind = argmax(z)
    print x[ind],y[ind]

    return x[ind],y[ind]

#Load in central influence function for RNR283 pixel 18
def RNR283(filename):
##    os.chdir('/Users/ryanallured/GoogleDrive/WFS/SystemAlignment/RNR283/140220')
    try:
        d = genfromtxt(filename)
    except:
        d = genfromtxt(filename,skip_header=5,delimiter='\t')   

    #Strip NaN rows/columns
    while sum(isnan(d[0]))==128:
        d = d[1:]
    while sum(isnan(d[-1]))==128:
        d = d[:-1]
    newsize = shape(d)[0]
    while sum(isnan(d[:,0]))==newsize:
        d = d[:,1:]
    while sum(isnan(d[:,-1]))==newsize:
        d = d[:,:-1]

    return -d

#Create x,y arrays to match RNR283, and fill height array
#using FEA model
def FEAmatch(filename):
    #Load in measurement and model
    d = RNR283(filename)*1000.
    xi,yi,zi = influence(24)

    #Create FEA model image array
    x,y = meshgrid(linspace(-50.,50.,shape(d)[1]),\
                   linspace(-50.,50.,shape(d)[0]))

    #Use griddata for interpolation, allows linear instead of nn
    npoints = (x,y)
    h = griddata((xi,yi),zi,npoints,method='linear')
    pdb.set_trace()

    #Set NaN region in model
    ind = isnan(d)
    h[ind] = NaN

    #Pad smaller dimension with NaN arrays
    s1,s2 = shape(d)
    pdb.set_trace()
    if s1<s2:
        for i in range(s2-s1):
            d = vstack((d,repeat(nan,s2)))
            h = vstack((h,repeat(nan,s2)))
    else:
        d = transpose(d)
        h = transpose(h)
        for i in range(s1-s2):
            d = vstack((d,repeat(nan,s1)))
            h = vstack((h,repeat(nan,s1)))
        d = transpose(d)
        h = transpose(h)

    #Scale peak to valley to influence function
    h = h - nanmin(h)
    h = h * (nanmax(d)-nanmin(d))/(nanmax(h))

    #Scale to match peaks
    h = h - nanmean(h) + nanmean(d)

    #Print RMS diff
    print sqrt(nanmean((h-d)**2))

    #Make plot showing overall results for central slice
    xarr = arange(shape(d)[0])
    ax1,ax2 = pltd(xarr,h[44],xarr,d[44]-h[44],\
            xlabel='Pixel',ylabel1='Height (nm)',ylabel2='Diff (nm)',\
                   title='RNR 283 Influence Measurement',label1='Model')
    ax1.plot(xarr,d[44],'b--',label='Measurement')
    ax1.legend(loc='lower center')
    

    return h,d

#Analyze RNR 283 Figure data for Pixel Pattern
def analyzepixelpattern(rnr283):
    #Load in RNR 283 figure wavefront
##    rnr283 = RNR283('/Users/ryanallured/GoogleDrive/WFS/SystemAlignment/RNR283/140625/140625FullAperture1.txt')
##    rnr283 = RNR283('/Users/ryanallured/GoogleDrive/WFS/SystemAlignment/CorningEagle/140714FullAperture1.txt')
##    rnr283 = RNR283('/Users/ryanallured/GoogleDrive/WFS/SystemAlignment/RNR283/140617/140618Repeatability.txt')
    #Pad smaller dimension with NaN arrays
    s1,s2 = shape(rnr283)
    if s1<s2:
        for i in range(s2-s1):
            rnr283 = vstack((rnr283,repeat(nan,s2)))
    else:
        rnr283 = transpose(rnr283)
        for i in range(s1-s2):
            rnr283 = vstack((rnr283,repeat(nan,s1)))
        rnr283 = transpose(rnr283)

    #Get center pixels of measurement and radius
    cx,cy,rad = locateimage(rnr283)
##    rad = 45.
    cx = round(cx)
    cy = round(cy)
    rad = floor(rad)

    #Fit Zernikes, remove first 90 from figure measurement, plus
    #N,N orders up to r=25
    r,m = zmodes(90)
    for i in range(13,26):
        r.append(i)
        r.append(i)
        m.append(i)
        m.append(-i)
    c = zcoeff(rnr283,cx=cx,cy=cy,rad=rad,r=r,m=m)
    resid = rnr283 - c[3]
    rem = c[3]

    pdb.set_trace()

    #Load in proper influence function
    f = '/Users/ryanallured/GoogleDrive/WFS/SystemAlignment/RNR283/ElectrodePattern/140708Summed90.txt'
    sinf = genfromtxt(f)

    #Filter Zernikes and remove
    cx,cy,rad = locateimage(sinf)
    c = zcoeff(sinf,cx=cx,cy=cy,rad=rad,r=r,m=m)
    sinf = sinf - c[3]

    return resid,sinf,rem

#Input residual ripple from figure measurement,
#make 2D PSD
#Reinterpolate on radial frequency grid
#Then integrate azimuthally and return 1D radial PSD
def radialPSD(rip,dx):
    #Pad rip input with 0s
    rip[isnan(rip)] = nanmean(rip)
    #Make 2D PSD
    freq = fft.fftshift(fft.fftfreq(shape(rip)[0],d=dx))
    psd2 = fft.fftshift(absolute(fft.fft2(rip))**2/size(freq)**4)
    #Divide by 2D frequency interval dfx*dfy
    psd2 = psd2/(abs(freq[0]-freq[1]))**2

    pdb.set_trace()

    #Convert to fx,fy,p arrays
    fx,fy = meshgrid(freq,freq)
    fx = fx.flatten()
    fy = fy.flatten()
    p = psd2.flatten()

    #Construct radial list of points
    nr = 1000.
    nt = 2000.
    frho = linspace(.01,.4,nr)
    ftheta = arange(0.,2*pi,2*pi/nt)
    frho,ftheta = meshgrid(frho,ftheta)
    frho = frho.flatten()
    ftheta = ftheta.flatten()
    fxp = frho*cos(ftheta)
    fyp = frho*sin(ftheta)

    #Run interpolation onto new grid
    points = transpose((fx,fy))
    newgrid = griddata(points,p,(fxp,fyp),method='linear')

    #Form 1D radial psd
    newgrid = reshape(newgrid,(nt,nr))
    frho = linspace(.01,.4,nr)
    psd1 = sum(newgrid,axis=0) * frho * (frho[1]-frho[0]) * (2*pi/nt)

    pdb.set_trace()

    return frho,psd1

#Cross correlation approach for RNR 283 data
def crosscorrelation():
    #0,0 corresponds to top left of wavefront looking at the sensor
    rnr283 = RNR283('/Users/ryanallured/GoogleDrive/WFS/SystemAlignment/RNR283/140625/140625FullAperture1.txt')
    #Load proper influence function
    s = genfromtxt('/Users/ryanallured/GoogleDrive/WFS/SystemAlignment/RNR283/ElectrodePattern/140804RNR283SummedInfRipple.txt')

    #Now you can do the cross correlation function - cell gap is at top of image
    snorm = nansum(s**2)
    rnrnorm = nansum(rnr283**2)
    
    
