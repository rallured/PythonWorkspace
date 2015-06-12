from numpy import *
from matplotlib.pyplot import *
import pdb
from utilities.plotting import pltd,mycontour
import glob

#Remove 2nd order poly
def remove2(l):
    l = l[invert(isnan(l))]
    x = arange(size(l))
    fit = polyfit(x,l,2)
    l = l - polyval(fit,x)

    return l

#Make PSD of an axial slice
#Dx just defines units of frequency
#Power is normalized to Parseval's theorem
#Divide by frequency interval when plotting
def axialPSD(l,dx,window=False,removeL2=True):
    l = l[invert(isnan(l))]
    length = dx * (size(l)-1)
    x = arange(size(l))*dx
    if removeL2 is True:
        fit = polyfit(x,l,2)
        l = l - polyval(fit,x)

    #Compute RMS
    rms = sqrt(mean(l**2))

    freq = fft.fftfreq(size(x),dx)
    if window is True:
        win = hanning(size(x))/sqrt(mean(hanning(size(x))**2))
    else:
        win = 1.
    spec = absolute(fft.fft(l*win))**2/size(l)**2
    rmsspec = sqrt(sum(spec))
    #Renormalize spectrum by Parseval's theorem
    if window is True:
        spec = spec * (rms/rmsspec)**2
    spec = 2*spec[freq>0]
    freq = freq[freq>0]

    return freq,spec

#Load in particular axial slice and take PSD
def analyzepsd(filename,n):
    #Load data
    try: 
        d = transpose(genfromtxt(filename,skip_header=5,delimiter='\t'))
    except:
        d = transpose(genfromtxt(filename,skip_header=6,delimiter='\t'))
        
    #Find longest slice, call that 25.4 mm, determine pixel size
    maxlength = 0
##    figure(1)
##    clf()
    for i in range(shape(d)[0]):
        l = d[i]
        l = l[invert(isnan(l))]
        maxlength = max([maxlength,size(l)])
##        plot(l)
    dx = 25.4 / (maxlength-1) #Pixel size in mm

    l = d[n]
    l = l[invert(isnan(l))]
    
    f,s = axialPSD(l,dx,window=False)

    return f,s

#Perform analysis on mid frequency wavefront
def analyzemidfreq(filename,plotpow=False,plotpsd=False,plotim=False,\
                   newfig=False,appendtable=False,\
                   dx=None,length=None,filtscrew=False):
    if type(filename) is str:
        #Load data
        try: 
            d = transpose(genfromtxt(filename))
        except:
            try:
                d = transpose(genfromtxt(filename,skip_header=5,delimiter='\t'))
            except:
                d = transpose(genfromtxt(filename,skip_header=6,delimiter='\t'))
    else:
        d = copy(filename)
        
    #Find longest slice, call that 25.4 mm, determine pixel size
    if length is not None:
        maxlength = 0
        for i in range(shape(d)[0]):
            l = d[i]
            l = l[invert(isnan(l))]
            maxlength= max([maxlength,size(l)])
    ##        plot(l)
        dx = length / (maxlength-1) #Pixel size in um

    #Make PSD for each axial slice with more than 50 valid pixels
    freq = []
    spec = []
    midpow = []
    midrange = []
    columns = []
    powarr = copy(d)

    sample = 'OP2S03'
    #Decompose filename to get sample and subaperture
##    sample = filename.split('Sample')[1][0]
##    if sample == '1':
##        sample = '12'
##        subaperture = filename.split('Sample')[1][2:].split('.')[0]
##    else:
##        subaperture = filename.split('Sample')[1][1:].split('.')[0]
##    if subaperture == 'BotLeft':
##        subaperture = '0'
##    elif subaperture == 'BotCenter':
##        subaperture = '1'
##    elif subaperture == 'BotRight':
##        subaperture = '2'
##    elif subaperture == 'MidLeft':
##        subaperture = '3'
##    elif subaperture == 'MidCenter':
##        subaperture = '4'
##    elif subaperture == 'MidRight':
##        subaperture = '5'
##    elif subaperture == 'TopLeft':
##        subaperture = '6'
##    elif subaperture == 'TopCenter':
##        subaperture = '7'
##    elif subaperture == 'TopRight':
##        subapertur = '8'

##    figure()
    if newfig is True and plotpsd is True:
        figure()
    for i in range(shape(d)[0]):
        l = d[i]
        l = l[invert(isnan(l))]
        if size(l) > 0:
            d[i][invert(isnan(d[i]))] = remove2(l)
            if filtscrew is True:
                d[i][invert(isnan(d[i]))] = \
                                filterscrew(d[i][invert(isnan(d[i]))])
        if size(l) > 89:
            columns.append(i)
            f,s = axialPSD(l,dx,window=True)
            if plotpsd is True:
                if i % 10 == 0:
                    loglog(f,s/(f[1]-f[0]),'r--',label=str(i))
            freq.append(f)
            spec.append(s)
            ind = logical_and(f>1/10.,f<1.)
            midpow.append(sqrt(sum(s[ind])))
            powarr[i,invert(isnan(d[i]))] = sqrt(sum(s[ind]))
            if appendtable is True:
                f = open('140616MayResults.txt','a')
                f.write(sample+'\t'+subaperture+'\t'+\
                        str(i)+'\t'+str(sqrt(sum(s[ind])))+'\n')
                f.close()
        else:
            powarr[i,invert(isnan(d[i]))] = nan


##    #Add repeatability PSD
##    if plotres is True:
##        #Get noise PSD
##        n = transpose(genfromtxt('/Users/ryanallured/GoogleDrive/WFS/SystemAlignment/Flat/140515Noise6.txt'\
##                                 ,skip_header=6,delimiter='\t'))
##        fn,sn = axialPSD(n[65],dx,window=True)
##        loglog(fn,sn/(fn[1]-fn[0]),label='Noise')
##        legend(loc='upper right')
##        title('Sample ' + sample + ' PSD')
##        xlabel('Frequency (1/mm)')
##        ylabel('Power ($\mu$m$^2$ mm)')
##        savefig(filename.split('.')[0]+'PSD.eps')
                                 

    #Make plots
    if plotpow is True:
        if newfig is True:
            figure()
        plot(array(columns),array(midpow)*1000.)
        xlabel('Azimuthal Pixel')
        ylabel('RMS Power (nm)')
##        savefig(filename.split('.')[0]+'MidPower.eps')
    if plotim is True:
        if newfig is True:
            figure()
##        mycontour(transpose(d),nlev=100)
        imshow(transpose(d*1000))
        colorbar()

    return freq,spec,d

#BLAH
def blahplot(fig=None):
    if fig is None:
        fig = figure()
    s = 1 
    for i in arange(9,0,-1):
        fig.add_subplot(3,3,s)
        analyzemidfreq(''+str(i)+'_2.txt',
                       plotpow=True,dx=25.4/114)
        ylim([0,35])
        s = s+1
    return fig

#May2014Measurements 9 sample plot
#Must be in MayMeasurements directory
def may2014plot(sample):
    fig = figure()
    f = glob.glob('*'+sample+'*.txt')
    s = 1
    for i in arange(9):
        fig.add_subplot(3,3,i+1)
        analyzemidfreq(f[i],plotim=True,dx=25.4/114)
    return fig

       
#Remove azimuthal sag for cylinder comparison
def removecyl(d,both=False):
    #Find fit parameters for central azimuthal slice
    sh = shape(d)
    cen = d[sh[1]/2]
    ind = invert(isnan(cen))
    x = arange(sh[0])
    fit = polyfit(x[ind],cen[ind],2)

    #Loop through azimuthal slices and remove quadratic
    for i in range(sh[0]):
        s = d[i]
        if sum(isnan(s))==sh[1]:
            continue
        #Identify valid range
        ind = invert(isnan(s))
        #Remove sag
        d[i] = s - polyval(fit,x)

    if both==True:
        #Find fit parameters for central axial slice
        cen = d[65]
        ind = invert(isnan(cen))
        x = arange(128.)
        fit = polyfit(x[ind],cen[ind],2)
        #Loop through axial slices and remove quadratic
        for i in range(128):
            s = d[i]
            if sum(isnan(s))==128:
                continue
            #Identify valid range
            ind = invert(isnan(s))
            #Remove sag
            d[i] = s - polyval(fit,x)

    return d

#Remove sag along the first axis
def removesag(d):
    for i in range(shape(d)[0]):
        d[i] = remove2(d[i])
    return d

#Filter frequencies higher than 1/3.175 mm^-1
def filterscrew(l):
    f = fft.fft(l)
    freq = fft.fftfreq(size(l),d=.05)
    f[where(abs(freq)>.25)]=0.
    return fft.ifft(f)

#Analyze section of NuStar slumped sample
def analyzeNuSection(d,axial1,axial2,az1,az2):
    d = copy(d[az1:az2,axial1:axial2])
    col,poww = analyzemidfreq(d,plotres=False,dx=.05,filtscrew=True)
    plot(col,poww)

#Analyze all measurements for a given sample
#Give strings defining the base filename for
#top, bottom, and middle measurements
#Plot the surfaces with cylinders removed
#Then overwrite the surface arrays with the midfrequency power
def analyzesample(botbase,midbase,topbase):
    #Construct array
    surfarr = zeros((128*3,128*3))
    powarr = copy(surfarr)
    
    surfarr[0:128,0:128],powarr[0:128,0:128] = analyzemidfreq(botbase+'Left.txt')
    surfarr[0:128,128:128*2],powarr[0:128,128:128*2] = analyzemidfreq(botbase+'Center.txt')
    surfarr[0:128,128*2:128*3],powarr[0:128,128*2:128*3] = analyzemidfreq(botbase+'Right.txt')

    surfarr[128:128*2,0:128],powarr[128:128*2,0:128] = analyzemidfreq(midbase+'Left.txt')
    surfarr[128:128*2,128:128*2],powarr[128:128*2,128:128*2] = analyzemidfreq(midbase+'Center.txt')
    surfarr[128:128*2,128*2:128*3],powarr[128:128*2,128*2:128*3] = analyzemidfreq(midbase+'Right.txt')

    surfarr[128*2:128*3,0:128],powarr[128*2:128*3,0:128] = analyzemidfreq(topbase+'Left.txt')
    surfarr[128*2:128*3,128:128*2],powarr[128*2:128*3,128:128*2] = analyzemidfreq(topbase+'Center.txt')
    surfarr[128*2:128*3,128*2:128*3],powarr[128*2:128*3,128*2:128*3] = analyzemidfreq(topbase+'Right.txt')

    powarr = powarr * 1000 #Convert to nm of RMS power

    return surfarr,powarr

