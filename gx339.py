import matplotlib.pyplot as plt
from plotting import *
import matplotlib
import numpy as np
import os
import scipy.optimize
import pdb
import glob
import pyfits
import scipy.stats
import sys

#Fe fixed to 1.592 based on fitting a constant
#nH fixed to 0.803 based on fitting a constant

##os.chdir('/Users/rallured/PythonWorkspace/GX339-4')
##data = np.genfromtxt('ResultsnHFeFixed3.txt',names='mflux,mfluxerr1,mfluxerr2 \
##,chisq,tin,tinerr1,tinerr2,dnorm,dnormerr1,dnormerr2,rin,rinerr1,rinerr2,phoind,phoinderr1,\
##phoinderr2,pnorm,pnormerr1,pnormerr2,xi,xierr1,xierr2,\
##rnorm,rnormerr1,rnormerr2')
##rsig2 = [19.6579,6.52981,24.9486,16.2301,60.1131,25,32.6835]
##rsig1 = [8.63143,3.3056,5.84274,6.65758,19.6691,16.0001,19.753]

##os.chdir('/Users/rallured/GX339-4/130324Results')
##i20e30 = np.genfromtxt('i20e30results.txt',names='chisq,nh,nherr1,nherr2,\
##tin,tinerr1,tinerr2,dnorm,dnormerr1,\
##dnormerr2,rin,rinerr1,rinerr2,phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
##,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
##i20e21 = np.genfromtxt('i20e21results.txt',names='chisq,nh,nherr1,nherr2,\
##tin,tinerr1,tinerr2,dnorm,dnormerr1,dnormerr2,rin,rinerr1,rinerr2,\
##phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
##,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
##i50e30 = np.genfromtxt('i50e30results.txt',names='chisq,nh,nherr1,nherr2,\
##tin,tinerr1,tinerr2,dnorm,dnormerr1,dnormerr2,rin,rinerr1,rinerr2,\
##phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
##,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
##i50e21 = np.genfromtxt('i50e21results.txt',names='chisq,nh,nherr1,nherr2,\
##tin,tinerr1,tinerr2,dnorm,dnormerr1,dnormerr2,rin,rinerr1,rinerr2,\
##phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
##,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
##i20e40 = np.genfromtxt('i20e40results.txt',names='chisq,nh,nherr1,nherr2,\
##tin,tinerr1,tinerr2,dnorm,dnormerr1,dnormerr2,rin,rinerr1,rinerr2,\
##phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
##,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
##i50e40 = np.genfromtxt('i50e40results.txt',names='chisq,nh,nherr1,nherr2,\
##tin,tinerr1,tinerr2,dnorm,dnormerr1,dnormerr2,rin,rinerr1,rinerr2,\
##phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
##,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
##os.chdir('..')

os.chdir('/Users/rallured/GX339-4/130409Results')
i20e30 = np.genfromtxt('i20e30Results.txt',names='\
tin,tinerr1,tinerr2,dnorm,dnormerr1,\
dnormerr2,rin,rinerr1,rinerr2,phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
i20e21 = np.genfromtxt('i20e21results.txt',names='\
tin,tinerr1,tinerr2,dnorm,dnormerr1,dnormerr2,rin,rinerr1,rinerr2,\
phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
i50e30 = np.genfromtxt('i50e30results2.txt',names='\
tin,tinerr1,tinerr2,dnorm,dnormerr1,dnormerr2,rin,rinerr1,rinerr2,\
phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
i50e21 = np.genfromtxt('i50e21results.txt',names='\
tin,tinerr1,tinerr2,dnorm,dnormerr1,dnormerr2,rin,rinerr1,rinerr2,\
phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
i20e40 = np.genfromtxt('i20e40results.txt',names='\
tin,tinerr1,tinerr2,dnorm,dnormerr1,dnormerr2,rin,rinerr1,rinerr2,\
phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
i50e40 = np.genfromtxt('i50e40results2.txt',names='\
tin,tinerr1,tinerr2,dnorm,dnormerr1,dnormerr2,rin,rinerr1,rinerr2,\
phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
os.chdir('..')

##data = np.genfromtxt('FinalResults2.txt',names='hflux,hfluxerr1,hfluxerr2, \
##sflux,sfluxerr1,sfluxerr2,chisq,tin,tinerr1,tinerr2,dnorm,dnormerr1,\
##dnormerr2,rin,rinerr1,rinerr2,phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
##,pnormerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2,tflux')

#Read in varied inclination/emissivity index data
##rin = np.zeros((8,10))
##rinerr1 = np.zeros((8,10))
##rinerr2 = np.zeros((8,10))
##chisq = np.zeros((8,10))
##i=0
##for file in glob.glob('spectrum7/loop??.out'):
##    temp = np.genfromtxt(file)
##    rin[i,:] = temp[0]
##    rinerr1[i,:] = temp[1]
##    rinerr2[i,:] = temp[2]
##    chisq[i,:] = temp[3]
##    i=i+1

#Read in PCA lightcurve data
pca = pyfits.open('lc/sourcelc_GX_339-4.lc')
pca = pca[1].data

#Observation times in MJD
#obstime = np.array([54245.,54261.,54285.,54299.,54733.,54897.,54908.,54908.5,55217.,55260.])
obstime = np.array([54245.,54261.,54285.,54299.,54897.,55217.,55260.])

##data = np.genfromtxt('FinalResults2.txt',names='sflux,hflux,tflux,\
##chisq,tin,tinerr1,tinerr2,dnorm,dnormerr1,dnormerr2,rin,rinerr1,rinerr2,\
##phoind,phoinderr1,phoinderr2,pnorm,pnormerr1,pnormerr2,xi,xierr1,xierr2,\
##rnorm,rnormerr1,rnormerr2')

#Read in data from end of May,beginning of June
##data = genfromtxt('FixedResults120606_2.txt',names='tin,tinerr1,tinerr2,dnorm,\
##    dnormerr1,dnormerr2,kind,kinderr1,kinderr2,rin,rinerr1,rinerr2,inc,incerr1,\
##    incerr2,pho,phoerr1,phoerr2,pnorm,pnormerr1,pnormerr2,fe,feerr1,feerr2,\
##    xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
##chisq = [1.2259,1.5398,0.783,1.1742,1.2013,1.5517,1.3253]

def addresults(arr,index,filename,chisq):
    arr2 = zeros((shape(arr)[0],shape(arr)[1]+1))
    arr2[0] = insert(arr[0],index,chisq)
    dadd = genfromtxt(filename)
    dadd = dadd.flatten()
    for i in arange(size(dadd))+1:
        arr2[i] = insert(arr[i],index,dadd[i-1])
    return arr2

def converterrors():
    data['mfluxerr1'] = abs(data['mflux'] - data['mfluxerr1'])
    data['mfluxerr2'] = abs(data['mflux'] - data['mfluxerr2'])
    data['tinerr1'] = abs(data['tin'] - data['tinerr1'])
    data['tinerr2'] = abs(data['tin'] - data['tinerr2'])
    data['dnormerr1'] = abs(data['dnorm'] - data['dnormerr1'])
    data['dnormerr2'] = abs(data['dnorm'] - data['dnormerr2'])
    data['rinerr1'] = abs(data['rin'] - data['rinerr1'])
    data['rinerr2'] = abs(data['rin'] - data['rinerr2'])
    data['phoinderr1'] = abs(data['phoind'] - data['phoinderr1'])
    data['phoinderr2'] = abs(data['phoind'] - data['phoinderr2'])
    data['pnormerr1'] = abs(data['pnorm'] - data['pnormerr1'])
    data['pnormerr2'] = abs(data['pnorm'] - data['pnormerr2'])
    data['xierr1'] = abs(data['xi'] - data['xierr1'])
    data['xierr2'] = abs(data['xi'] - data['xierr2'])
    data['rnormerr1'] = abs(data['rnorm'] - data['rnormerr1'])
    data['rnormerr2'] = abs(data['rnorm'] - data['rnormerr2'])

def chooseasymerrors(const,pts,lowerror,higherror):
    asymerror = np.zeros(np.size(pts))
    for i in range(np.size(pts)):
        if pts[i] < const:
            asymerror[i] = higherror[i]
        else:
            asymerror[i] = lowerror[i]
    #asymerror = np.abs(asymerror-pts)
    return asymerror

def fconst(x,c):
    return np.repeat(c,np.size(x))

#Perrs must be relative (+- values)
def fitconst(p,perr1,perr2,carr):

    chi = np.zeros(np.size(carr))
    for i in range(np.size(carr)):
        err = chooseasymerrors(carr[i],p,perr1,perr2)
        chi[i] = sum((p-carr[i])**2/(err)**2)/(np.size(p)-1)
        
    return chi

def makeplots(d):
    plt.ion()
    plt.hold(True)
##    matplotlib.rcParams['font.weight']='bolder'
##    matplotlib.rc('text',usetex=True)

    x = np.arange(1,8)

    os.chdir('/Users/rallured/GX339-4/130324Results/i50e21/')

##    #Make flux plot
##    plt.clf()
##    #plt.errorbar(x,data['mflux'],yerr=[data['mfluxerr1'],data['mfluxerr2']]\
##    #             ,fmt='.',ms=12)
##    plt.plot(x,data['tflux'],'r.')
####    for i in range(7):
####        plt.plot([x[i],x[i]],[data['mfluxerr2'][i],\
####                              data['mfluxerr1'][i]],'r')
##    plt.title('Unabsorbed Flux vs. Spectrum Number')
##    plt.xlabel('Spectrum Number')
##    plt.ylabel('Flux (photons/cm$^2$/s)')
##    plt.xlim([0,8])
##    plt.savefig('Flux.eps')

    #Make chisq plot
    plt.clf()
    plt.plot(x,d['chisq'],'.',ms=12)
    plt.title('Chi Squared vs. Spectrum Number')
    plt.xlabel('Spectrum Number')
    plt.ylabel('Reduced Chi Squared')
    plt.xlim([0,11])
    plt.savefig('Chisq.eps')

    #Make nH plot
    plt.clf()
    plt.errorbar(x,d['nh'],yerr=[d['nherr1'],d['nherr2']]\
                 ,fmt='.',ms=12)
    plt.title('Hydrogen Column Density vs. Spectrum Number')
    plt.xlabel('Spectrum Number')
    plt.ylabel('Hydrogen Column Density (atoms/cm$^2$)')
    plt.xlim([0,11])
    plt.savefig('nH.eps')

    #Make Tin plot
    plt.clf()
    ind = [0,1,2,3,4,6]
    print ind
    plt.plot(x[ind],d['tin'][ind],'r.')
    for i in ind:
        plt.plot([x[i],x[i]],[d['tinerr2'][i],\
                              d['tinerr1'][i]],'r')
    plt.title('Inner Temperature vs. Spectrum Number')
    plt.xlabel('Spectrum Number')
    plt.ylabel('Inner Temperature (keV)')
    plt.xlim([0,11])
    plt.savefig('tin.eps')

    #Make dnorm plot
    plt.clf()
    plt.semilogy(x[ind],d['dnorm'][ind],'r.')
    for i in ind:
        plt.plot([x[i],x[i]],[d['dnormerr2'][i],\
                              d['dnormerr1'][i]],'r')
    #plt.errorbar(x,data['dnorm'],yerr=[data['dnormerr1']\
    #              ,data['dnormerr2']],fmt='.',ms=12)
    plt.title('Disk Normalization vs. Spectrum Number')
    plt.xlabel('Spectrum Number')
    plt.ylabel('Disk Normalization')
    plt.xlim([0,11])
    plt.savefig('dnorm.eps')

    #Make kind plot
##    plt.clf()
##    plt.plot(x,data['kind'],'r.')
##    for i in range(7):
##        plt.plot([x[i],x[i]],[data['kinderr2'][i],\
##                              data['kinderr1'][i]],'r')
##    plt.title('Kdblur Emissivity Index vs. Spectrum Number')
##    plt.xlabel('Spectrum Number')
##    plt.ylabel('Emissivity Index')
##    plt.xlim([0,8])
##    plt.ylim([1.9,3.1])
##    plt.savefig('kind.eps')

    #Make inc plot
##    plt.clf()
##    plt.plot(x,data['inc'],'r.')
##    for i in range(7):
##        plt.plot([x[i],x[i]],[data['incerr2'][i],\
##                              data['incerr1'][i]],'r')
##    plt.title('Disk Inclination vs. Spectrum Number')
##    plt.xlabel('Spectrum Number')
##    plt.ylabel('Disk Inclination (deg)')
##    plt.xlim([0,8])
##    plt.savefig('inc.eps')

    #Make rin plot
    plt.clf()
    plt.plot(x,d['rin'],'r.')
    for i in range(7):
        plt.plot([x[i],x[i]],[d['rinerr2'][i],\
                              d['rinerr1'][i]],'r')
    plt.title('Inner Radius vs. Spectrum Number')
    plt.xlabel('Spectrum Number')
    plt.ylabel('Kdblur Inner Radius ($GM/c^2$)')
    plt.xlim([0,11])
    plt.savefig('rin.eps')

    #Make Phoind plot
    plt.clf()
    #plt.errorbar(x,data['phoind'],yerr=[data['phoinderr1'],data['phoinderr2']]\
    #             ,fmt='.',ms=12)
    plt.plot(x,d['phoind'],'r.')
    for i in range(7):
        plt.plot([x[i],x[i]],[d['phoinderr2'][i],\
                              d['phoinderr1'][i]],'r')
    plt.title('Photon Index vs. Spectrum Number')
    plt.xlabel('Spectrum Number')
    plt.ylabel('Photon Index')
    plt.xlim([0,11])
    plt.savefig('phoind.eps')

    #Make Pnorm plot
    plt.clf()
    #plt.errorbar(x,data['pnorm'],yerr=[data['pnormerr1'],data['pnormerr2']]\
    #             ,fmt='.',ms=12)
    plt.plot(x,d['pnorm'],'r.')
    for i in range(7):
        plt.plot([x[i],x[i]],[d['pnormerr2'][i],\
                              d['pnormerr1'][i]],'r')
    plt.title('Power Law Normalization vs. Spectrum Number')
    plt.xlabel('Spectrum Number')
    plt.ylabel('Power Law Normalization')
    plt.xlim([0,11])
    plt.savefig('pnorm.eps')

    #Make fe plot
    plt.clf()
    plt.plot(x,d['fe'],'r.')
    for i in range(7):
        plt.plot([x[i],x[i]],[d['feerr2'][i],\
                              d['feerr1'][i]],'r')
    plt.title('Iron Abundance vs. Spectrum Number')
    plt.xlabel('Spectrum Number')
    plt.ylabel('Iron Abundance (solar)')
    plt.xlim([0,11])
    plt.savefig('iron.eps')

    #Make xi plot
    plt.clf()
    plt.plot(x,d['xi'],'r.')
    for i in range(7):
        plt.plot([x[i],x[i]],[d['xierr2'][i],\
                              d['xierr1'][i]],'r')
    plt.title('Ionization Parameter vs. Spectrum Number')
    plt.xlabel('Spectrum Number')
    plt.ylabel('Ionization Parameter (erg cm/s)')
    plt.xlim([0,11])
    plt.savefig('Xi.eps')

    #Make reflion norm plot
    plt.clf()
    #plt.errorbar(x,data['rnorm'],yerr=[data['rnormerr1'],data['rnormerr2']]\
    #             ,fmt='.',ms=12)
    plt.plot(x,d['rnorm'],'r.')
    for i in range(7):
        plt.plot([x[i],x[i]],[d['rnormerr2'][i],\
                              d['rnormerr1'][i]],'r')
    plt.title('Reflionx Normalization vs. Spectrum Number')
    plt.xlabel('Spectrum Number')
    plt.ylabel('Reflionx Normalization')
    plt.xlim([0,11])
    plt.savefig('Rnorm.eps')

    #Make varied inclination/emissivity index plot
##    plt.clf()
##    fig = figure(figsize=(6,18))
##    emis = np.arange(3.0,2.0,-.1)
##    incl = np.arange(10,90,10)
##    colr = ['r','b','g','y','c','m','k','chartreuse']
##    for i in range(8):
##        a = fig.add_subplot(8,1,i+1)
##        if i!=7:
##            a.set_xticklabels('',visible=False)
##        else:
##            plt.xlabel('Emissivity Index')
##        a.plot(emis,rin[i,:],'.',color=colr[i],\
##                     label='i='+str(incl[i]))
##        plt.title('Inclination = ' + str(incl[i]))
##        plt.ylabel('$R_{in}$ ($GM/c^2$)')
##        for m in range(size(emis)):
##            plt.plot([emis[m],emis[m]],[rinerr1[i,m],rinerr2[i,m]],\
##                     '-',color=colr[i])
##    #plt.xlim([2.0,3.4])
##    #plt.legend(loc='upper right')
##    #plt.title('$R_{in}$ vs. Emissivity Index')
##    #plt.xlabel('Emissivity Index')
##    plt.savefig('InclEmisSubplots.eps')
##    plt.clf()
##    mycontour(chisq,1,x=emis,y=incl)
##    plt.title('$\chi^2$ vs. Emissivity Index and Inclination')
##    plt.xlabel('Emissivity Index')
##    plt.ylabel('Inclination (deg)')
##    savefig('InclEmisChisq.eps')

    #Make Rin and lightcurve panel plot
    plt.clf()
    fig = figure(1)
    ax = plt.subplot(2,1,1)
    fig.subplots_adjust(hspace=.12)
    plt.setp(ax.get_xticklabels(),visible=False)
    plt.plot(obstime-50000,d['rin'],'r.')
    for i in range(7):
        plt.plot(np.array([obstime[i],obstime[i]])-50000,\
                 [d['rinerr2'][i],d['rinerr1'][i]],'r')
    plt.arrow(4733,65,0,10,head_length=5,head_width=20,\
              color='green',linewidth=1.5)
    plt.plot([4733-20,4733+20],[65,65],'g-')
    plt.title('$R_{in}$ and PCA Rate vs. Observation Time')
    plt.ylabel('$R_{in}$ ($GM/c^2$)')
    plt.xlim([4100,5400])
    plt.ylim([0,80])
    plt.subplot(2,1,2)
    plt.errorbar(pca['time']-50000,pca['rate'],pca['error'],fmt='.')
    plt.xlabel('Observation Time (MJD-50000)')
    plt.ylabel('PCA Rate (cts/sec)')
    for i in range(7):
        plt.plot(np.array([obstime[i],obstime[i]])-50000,\
                 [-1000,7000],'k--',linewidth=1)
    plt.plot([4733,4733],[-1000,7000],'k--',linewidth=1)
    plt.xlim([4100,5400])
    plt.savefig('RinLC.eps')

    #Make Rin Dnorm comparison
    n_dnorm = sqrt(d['dnorm'])/sqrt(min(d['dnorm'][ind]))
    n_dnormerr1 = sqrt(d['dnormerr1'])/sqrt(min(d['dnorm'][ind]))
    n_dnormerr2 = sqrt(d['dnormerr2'])/sqrt(min(d['dnorm'][ind]))
    n_rin = d['rin']/min(d['rin'])
    n_rinerr1 = d['rinerr1']/min(d['rin'])
    n_rinerr2 = d['rinerr2']/min(d['rin'])
    plt.clf()
##    plt.plot(x[ind],n_dnorm[ind],'b.',label='DiskNorm')
##    plt.plot(x,n_rin,'r.',label='R_in')
##    #Plot out sqrt(dnorm) error bars
##    for i in ind:
##        plt.plot([x[i],x[i]],[n_dnormerr1[i],n_dnormerr2[i]],'b')
##    for i in range(7):
##        plt.plot([x[i],x[i]],[n_rinerr1[i],n_rinerr2[i]],'r')
##    plt.xlim([0,8])
    plt.errorbar(x[ind],n_dnorm[ind],yerr=[n_dnorm[ind]-n_dnormerr1[ind],\
        n_dnormerr2[ind]-n_dnorm[ind]],label='D$_{norm}$',fmt='.',capsize=10)
    plt.errorbar(x,n_rin,yerr=[n_rin-n_rinerr1,n_rinerr2-n_rin],\
                 label='R$_{in}$',fmt='.',capsize=10)
    plt.xlim([0,11])
    plt.legend(loc='upper right')
    plt.title('Inner Radius Comparison')
    plt.xlabel('Spectrum Number')
    plt.ylabel('Inner Radius (relative to spectrum 1)')
    plt.savefig('MethodComparison.eps')

    #Rin Dnorm Scatterplot
    plt.clf()
    plt.errorbar(n_dnorm,n_rin,xerr=[n_dnormerr1,n_dnormerr2],\
                 yerr=[n_rinerr1,n_rinerr2],fmt='b.',ecolor='r')
    plt.xlabel('Sqrt(D$_{norm}$)')
    plt.ylabel('R$_{in}$')
    plt.title('Inner Radius Comparison')
    plt.savefig('RinvsDnorm.eps')

##    #Make luminosity-radius correlation plots
##    timelag = 100.
##    lum = np.zeros(7)
##    lum2 = np.zeros(7)
##    for i in range(7):
##        diff = pca['time']-(obstime[i]-timelag)
##        close = where(abs(diff)==min(abs(diff)))[0]
##        if pca['time'][close] > obstime[i]-timelag:
##            close = close-1
##        slope = (pca['rate'][close+1]-pca['rate'][close])\
##                /(pca['time'][close+1]-pca['time'][close])
##        lum[i] = pca['rate'][close] + slope*\
##                 (obstime[i]-timelag-pca['time'][close])
##        diff = pca['time']-(obstime[i])
##        close = where(abs(diff)==min(abs(diff)))[0]
##        if pca['time'][close] > obstime[i]:
##            close = close-1
##        slope = (pca['rate'][close+1]-pca['rate'][close])\
##                /(pca['time'][close+1]-pca['time'][close])
##        lum2[i] = pca['rate'][close] + slope*\
##                 (obstime[i]-pca['time'][close])
##    print lum
##    plt.clf()
##    fig = figure(1)
##    ax = plt.subplot(2,1,1)
##    fig.subplots_adjust(hspace=.12)
##    plt.setp(ax.get_xticklabels(),visible=False)
##    plt.semilogx(lum,data['rin'],'b.')
##    plt.errorbar(lum,data['rin'],[data['rin']-data['rinerr1'],\
##                                  data['rinerr2']-data['rin']],fmt='b.')
##    plt.title('Inner Radius vs. PCA Rate')
##    plt.ylabel('Inner Radius $(GM/c^2)$')
##    plt.subplot(2,1,2)
##    plt.semilogx(lum2,data['rin'],'g.')
##    plt.errorbar(lum2,data['rin'],[data['rin']-data['rinerr1'],\
##                                  data['rinerr2']-data['rin']],fmt='g.')
##    plt.xlabel('PCA Rate (cts/sec)')
##    plt.ylabel('Inner Radius $(GM/c^2)$')
##    plt.savefig('TimeLag.eps')

##feerr1 = np.array([1.18858,1.24899,0.861168,1.66394,0.918209,1.0104,1.52184])
##feerr2 = np.array([1.71914,2.35114,1.46451,2.71313,1.2313,1.72816,1.99669])
##fe = np.array([1.45063,1.85211,1.05759,2.0012,1.05211,1.37111,1.80326]) #femin = 1.391
##
##nherr1 = np.array([.766231,.854993,0.882055,0.755188,0.738646,0.676568,0.817005])
##nherr2 = np.array([.861275,1.02513,1.12428,0.924733,0.833834,0.730035,0.914126])
##nh = np.array([.809033,.926953,1.0107,0.865119,0.788498,0.699483,0.861654]) #nhmin = .788

##nherr1 = np.array([.738665,.844859,.755441,.749304,.710284,.679441,.795462])
##nherr2 = np.array([.825935,.984232,1.03361,.906814,.824892,.725362,.873627])
##nh = np.array([.789327,.921819,.901992,.823607,.776283,.701128,.828815])
##
##feerr1 = np.array([1.59024,1.48572,1.47554,1.88317,.943167,.83193,1.40815])
##feerr2 = np.array([2.52451,3.28566,4.66327,2.51711,1.44039,1.35983,1.87274])
##fe = np.array([1.99854,2.00058,2.33495,1.99699,1.09104,1.06684,1.64529])
#Yields an nH of .766 and a Fe of 1.804

#Normalize rin and dnorm for comparison
def rindnorm(d):
    ind = [0,1,2,3,6]
    n_dnorm = sqrt(d['dnorm'])#/sqrt(min(d['dnorm'][ind]))
    n_dnormerr1 = sqrt(d['dnormerr1'])#/sqrt(min(d['dnorm'][ind]))
    n_dnormerr2 = sqrt(d['dnormerr2'])#/sqrt(min(d['dnorm'][ind]))
    n_rin = d['rin']#/min(d['rin'])
    n_rinerr1 = d['rinerr1']#/min(d['rin'])
    n_rinerr2 = d['rinerr2']#/min(d['rin'])
    return n_dnorm, n_dnormerr1, n_dnormerr2, n_rin, n_rinerr1, n_rinerr2

#Make sqrt(dnorm) rin comparison for all 4 rounds
def dnormrincomp():
    fig = plt.figure(1)
    plt.clf()

    matplotlib.rcParams['axes.titlesize']=12
    matplotlib.rcParams['axes.labelsize']=12
    
    n_dnorm, n_dnormerr1, n_dnormerr2, n_rin, n_rinerr1, n_rinerr2 = rindnorm(i20e30)
    fig.add_subplot(3,2,1)
    plt.loglog(n_dnorm,n_rin,visible=False)
    plt.errorbar(n_dnorm,n_rin,xerr=[n_dnorm-n_dnormerr1,n_dnormerr2-n_dnorm],\
                 yerr=[n_rin-n_rinerr1,n_rinerr2-n_rin],fmt='b.',ecolor='r')
    plt.xlabel('i=20, e=3.0 $\sqrt{D_{norm}}$')
    plt.ylabel('i=20, e=3.0 R$_{in}$')
##    plt.title('Inner Radius Comparison')

    n_dnorm, n_dnormerr1, n_dnormerr2, n_rin, n_rinerr1, n_rinerr2 = rindnorm(i20e21)
    fig.add_subplot(3,2,2)
    plt.loglog(n_dnorm,n_rin,visible=False)
    plt.errorbar(n_dnorm,n_rin,xerr=[n_dnorm-n_dnormerr1,n_dnormerr2-n_dnorm],\
                 yerr=[n_rin-n_rinerr1,n_rinerr2-n_rin],fmt='b.',ecolor='r')
    plt.xlabel('i=20, e=2.1 $\sqrt{D_{norm}}$')
    plt.ylabel('i=20, e=2.1 R$_{in}$')

    n_dnorm, n_dnormerr1, n_dnormerr2, n_rin, n_rinerr1, n_rinerr2 = rindnorm(i50e30)
    fig.add_subplot(3,2,3)
    plt.loglog(n_dnorm,n_rin,visible=False)
    plt.errorbar(n_dnorm,n_rin,xerr=[n_dnorm-n_dnormerr1,n_dnormerr2-n_dnorm],\
                 yerr=[n_rin-n_rinerr1,n_rinerr2-n_rin],fmt='b.',ecolor='r')
    plt.xlabel('i=50, e=3.0 $\sqrt{D_{norm}}$')
    plt.ylabel('i=50, e=3.0 R$_{in}$')

    n_dnorm, n_dnormerr1, n_dnormerr2, n_rin, n_rinerr1, n_rinerr2 = rindnorm(i50e21)
    fig.add_subplot(3,2,4)
    plt.loglog(n_dnorm,n_rin,visible=False)
    plt.errorbar(n_dnorm,n_rin,xerr=[n_dnorm-n_dnormerr1,n_dnormerr2-n_dnorm],\
                 yerr=[n_rin-n_rinerr1,n_rinerr2-n_rin],fmt='b.',ecolor='r')
    plt.xlabel('i=50, e=2.1 $\sqrt{D_{norm}}$')
    plt.ylabel('i=50, e=2.1 R$_{in}$')

    n_dnorm, n_dnormerr1, n_dnormerr2, n_rin, n_rinerr1, n_rinerr2 = rindnorm(i50e40)
    fig.add_subplot(3,2,5)
    plt.loglog(n_dnorm,n_rin,visible=False)
    plt.errorbar(n_dnorm,n_rin,xerr=[n_dnorm-n_dnormerr1,n_dnormerr2-n_dnorm],\
                 yerr=[n_rin-n_rinerr1,n_rinerr2-n_rin],fmt='b.',ecolor='r')
    plt.xlabel('i=50, e=4.0 $\sqrt{D_{norm}}$')
    plt.ylabel('i=50, e=4.0 R$_{in}$')

    n_dnorm, n_dnormerr1, n_dnormerr2, n_rin, n_rinerr1, n_rinerr2 = rindnorm(i20e40)
    fig.add_subplot(3,2,6)
    plt.loglog(n_dnorm,n_rin,visible=False)
    plt.errorbar(n_dnorm,n_rin,xerr=[n_dnorm-n_dnormerr1,n_dnormerr2-n_dnorm],\
                 yerr=[n_rin-n_rinerr1,n_rinerr2-n_rin],fmt='b.',ecolor='r')
    plt.xlabel('i=20, e=4.0 $\sqrt{D_{norm}}$')
    plt.ylabel('i=20, e=4.0 R$_{in}$')

    fig.suptitle('R$_{in}$ vs. $\sqrt{D_{norm}}$',size=18)

    matplotlib.rcParams['axes.titlesize']=18
    matplotlib.rcParams['axes.labelsize']=16

#Take two rounds and normalize their inner radii
def normrin(d1,d2):
##    x = array([0,1,2,3,4,5,6])
##    norm = min(d1['rin'][x])
##    n1_rin = d1['rin'][x]/norm
##    n1_rinerr1 = n1_rin-d1['rinerr1'][x]/norm
##    n1_rinerr2 = d1['rinerr2'][x]/norm-n1_rin
##    norm = min(d2['rin'][x])
##    n2_rin = d2['rin'][x]/norm
##    n2_rinerr1 = n2_rin-d2['rinerr1'][x]/norm
##    n2_rinerr2 = d2['rinerr2'][x]/norm-n2_rin

    x = array([0,1,2,3,4,5])
    norm = 1.#min(d1['dnorm'][x])
    n1_rin = d1['dnorm'][x]/norm
    n1_rinerr1 = n1_rin-d1['dnormerr1'][x]/norm
    n1_rinerr2 = d1['dnormerr2'][x]/norm-n1_rin
    norm = 1.#min(d2['dnorm'][x])
    n2_rin = d2['dnorm'][x]/norm
    n2_rinerr1 = n2_rin-d2['dnormerr1'][x]/norm
    n2_rinerr2 = d2['dnormerr2'][x]/norm-n2_rin
    return n1_rin, n1_rinerr1, n1_rinerr2, n2_rin, n2_rinerr1, n2_rinerr2

#Make rin comparison plots between different emissivity indices
#at same inclination
def ecomp():
    matplotlib.rcParams['axes.titlesize']=12
    matplotlib.rcParams['axes.labelsize']=12
    
    n1_rin, n1_rinerr1, n1_rinerr2, n2_rin, n2_rinerr1, n2_rinerr2\
            = normrin(i20e30,i50e30)
    fig = plt.figure(1)
    plt.clf()
    fig.add_subplot(421)
    plt.loglog(n1_rin,n2_rin,visible=False)
    plt.errorbar(n1_rin,n2_rin,xerr=[n1_rinerr1,n1_rinerr2],\
                yerr=[n2_rinerr1,n2_rinerr2],ecolor='r',fmt='b.')
##    plt.title('Relative R$_{in}$ Change Comparison: i var, e=3.0')
    plt.xlabel('i=20, e=3.0 $\sqrt{D_{norm}}$')
    plt.ylabel('i=50, e=3.0 $\sqrt{D_{norm}}$')
##    plt.xlabel('i=20, e=3.0 R$_{in}$ (arbitrary units)')
##    plt.ylabel('i=50, e=3.0 R$_{in}$ (arbitrary units)')

    n1_rin, n1_rinerr1, n1_rinerr2, n2_rin, n2_rinerr1, n2_rinerr2\
            = normrin(i20e21,i50e21)
    fig.add_subplot(422)
    plt.loglog(n1_rin,n2_rin,visible=False)
    plt.errorbar(n1_rin,n2_rin,xerr=[n1_rinerr1,n1_rinerr2],\
                yerr=[n2_rinerr1,n2_rinerr2],ecolor='r',fmt='b.')
##    plt.title('Relative R$_{in}$ Change Comparison: i var, e=2.1')
    plt.xlabel('i=20, e=2.1 $\sqrt{D_{norm}}$')
    plt.ylabel('i=50, e=2.1 $\sqrt{D_{norm}}$')
##    plt.xlabel('i=20, e=2.1 R$_{in}$ (arbitrary units)')
##    plt.ylabel('i=50, e=2.1 R$_{in}$ (arbitrary units)')

    n1_rin, n1_rinerr1, n1_rinerr2, n2_rin, n2_rinerr1, n2_rinerr2\
            = normrin(i50e21,i50e30)
    fig.add_subplot(423)
    plt.loglog(n1_rin,n2_rin,visible=False)
    plt.errorbar(n1_rin,n2_rin,xerr=[n1_rinerr1,n1_rinerr2],\
                yerr=[n2_rinerr1,n2_rinerr2],ecolor='r',fmt='b.')
##    plt.title('Relative R$_{in}$ Change Comparison: i=50, e var')
    plt.xlabel('i=50, e=2.1 $\sqrt{D_{norm}}$')
    plt.ylabel('i=50, e=3.0 $\sqrt{D_{norm}}$')
##    plt.xlabel('i=50, e=2.1 R$_{in}$ (arbitrary units)')
##    plt.ylabel('i=50, e=3.0 R$_{in}$ (arbitrary units)')
    
    n1_rin, n1_rinerr1, n1_rinerr2, n2_rin, n2_rinerr1, n2_rinerr2\
            = normrin(i20e21,i20e30)
    fig.add_subplot(424)
    plt.loglog(n1_rin,n2_rin,visible=False)
    plt.errorbar(n1_rin,n2_rin,xerr=[n1_rinerr1,n1_rinerr2],\
                yerr=[n2_rinerr1,n2_rinerr2],ecolor='r',fmt='b.')
##    plt.title('Relative R$_{in}$ Change Comparison: i=20, e var')
    plt.xlabel('i=20, e=2.1 $\sqrt{D_{norm}}$')
    plt.ylabel('i=20, e=3.0 $\sqrt{D_{norm}}$')
##    plt.xlabel('i=20, e=2.1 R$_{in}$ (arbitrary units)')
##    plt.ylabel('i=20, e=3.0 R$_{in}$ (arbitrary units)')

    n1_rin, n1_rinerr1, n1_rinerr2, n2_rin, n2_rinerr1, n2_rinerr2\
            = normrin(i20e21,i50e30)
    fig.add_subplot(425)
    plt.loglog(n1_rin,n2_rin,visible=False)
    plt.errorbar(n1_rin,n2_rin,xerr=[n1_rinerr1,n1_rinerr2],\
                yerr=[n2_rinerr1,n2_rinerr2],ecolor='r',fmt='b.')
##    plt.title('Relative R$_{in}$ Change Comparison: i var, e var')
    plt.xlabel('i=20, e=2.1 $\sqrt{D_{norm}}$')
    plt.ylabel('i=50, e=3.0 $\sqrt{D_{norm}}$')
##    plt.xlabel('i=20, e=2.1 R$_{in}$ (arbitrary units)')
##    plt.ylabel('i=50, e=3.0 R$_{in}$ (arbitrary units)')

    n1_rin, n1_rinerr1, n1_rinerr2, n2_rin, n2_rinerr1, n2_rinerr2\
            = normrin(i20e30,i50e21)
    fig.add_subplot(426)
    plt.loglog(n1_rin,n2_rin,visible=False)
    plt.errorbar(n1_rin,n2_rin,xerr=[n1_rinerr1,n1_rinerr2],\
                yerr=[n2_rinerr1,n2_rinerr2],ecolor='r',fmt='b.')
##    plt.title('Relative R$_{in}$ Change Comparison: i var, e var')
    plt.xlabel('i=20, e=3.0 $\sqrt{D_{norm}}$')
    plt.ylabel('i=50, e=2.1 $\sqrt{D_{norm}}$')
##    plt.xlabel('i=20, e=3.0 R$_{in}$ (arbitrary units)')
##    plt.ylabel('i=50, e=2.1 R$_{in}$ (arbitrary units)')

    n1_rin, n1_rinerr1, n1_rinerr2, n2_rin, n2_rinerr1, n2_rinerr2\
            = normrin(i20e30,i20e40)
    fig.add_subplot(427)
    plt.loglog(n1_rin,n2_rin,visible=False)
    plt.errorbar(n1_rin,n2_rin,xerr=[n1_rinerr1,n1_rinerr2],\
                yerr=[n2_rinerr1,n2_rinerr2],ecolor='r',fmt='b.')
##    plt.title('Relative R$_{in}$ Change Comparison: i var, e var')
    plt.xlabel('i=20, e=3.0 $\sqrt{D_{norm}}$')
    plt.ylabel('i=20, e=4.0 $\sqrt{D_{norm}}$')
##    plt.xlabel('i=20, e=3.0 R$_{in}$ (arbitrary units)')
##    plt.ylabel('i=50, e=2.1 R$_{in}$ (arbitrary units)')

    n1_rin, n1_rinerr1, n1_rinerr2, n2_rin, n2_rinerr1, n2_rinerr2\
            = normrin(i50e30,i50e40)
    fig.add_subplot(428)
    plt.loglog(n1_rin,n2_rin,visible=False)
    plt.errorbar(n1_rin,n2_rin,xerr=[n1_rinerr1,n1_rinerr2],\
                yerr=[n2_rinerr1,n2_rinerr2],ecolor='r',fmt='b.')
##    plt.title('Relative R$_{in}$ Change Comparison: i var, e var')
    plt.xlabel('i=50, e=3.0 $\sqrt{D_{norm}}$')
    plt.ylabel('i=50, e=4.0 $\sqrt{D_{norm}}$')
##    plt.xlabel('i=20, e=3.0 R$_{in}$ (arbitrary units)')
##    plt.ylabel('i=50, e=2.1 R$_{in}$ (arbitrary units)')

    fig.suptitle('Relative $\sqrt{D_{norm}}$ Change Comparison',size=18)
##    fig.suptitle('Relative R$_{in}$ Change Comparison',size=18)

    matplotlib.rcParams['axes.titlesize']=18
    matplotlib.rcParams['axes.labelsize']=16


def combinerin():
    lowr = np.zeros(7)
    highr = copy(lowr)
    lowr2 = copy(lowr)
    highr2 = copy(lowr)
    for i in range(7):
        lowr[i] = np.min([i20e30['rinerr1'][i],i20e21['rinerr1'][i],\
                          i50e30['rinerr1'][i],i50e21['rinerr1'][i]])
        highr[i] = np.max([i20e30['rinerr2'][i],i20e21['rinerr2'][i],\
                          i50e30['rinerr2'][i],i50e21['rinerr2'][i]])
        lowr2[i] = np.min([i50e30['rinerr1'][i],i50e21['rinerr1'][i]])
        highr2[i] = np.max([i50e30['rinerr2'][i],i50e21['rinerr2'][i]])
    plt.clf()
    plt.semilogy(arange(1,8),repeat(50,7),visible=False)    
    plt.errorbar(arange(1,8),repeat(50,7),\
                 yerr=[50-lowr,highr-50],visible=False)
    plt.errorbar(arange(1,8),repeat(50,7),\
                 yerr=[50-lowr2,highr2-50],visible=False)
    plt.title('Combined Inner Radius Errors')
    plt.xlabel('Spectrum \#')
    plt.ylabel('Inner Radius ($GM/c^2$)')
    plt.xlim([0,8])
    print lowr
    print highr

#Make Rin vs Time plot
def rinlc():
    plt.clf()
    fig = figure(1)
    ax1 = plt.subplot(5,1,1)
    ax1.set_xticklabels('',visible=False)
    plt.title('$R_{in}$ vs. Observation Time (Kdblur)')
    x = [0,1,2,3,4,6]
    ax1.errorbar(obstime[x]-50000,i20e30['rin'][x],yerr=[i20e30['rin'][x]-\
        i20e30['rinerr1'][x],i20e30['rinerr2'][x]-i20e30['rin'][x]],fmt='r.')
    ax1.errorbar(obstime[5]-50000,mean([i20e30['rinerr1'][5],100]),\
                 yerr=100-mean([i20e30['rinerr1'][5],100]),\
                 uplims=True,fmt='r')
    ax1.set_xlim([4100,5400])
    ax1.set_ylim([0,100])
    ax1.set_ylabel('i=20,e=3.0')
    ax1.set_yticklabels(('0','20','40','60','80','100'))
    
    ax2 = plt.subplot(5,1,2)
    ax2.set_xticklabels('',visible=False)
    ax2.errorbar(obstime[x]-50000,i20e21['rin'][x],yerr=[i20e21['rin'][x]-\
        i20e21['rinerr1'][x],i20e21['rinerr2'][x]-i20e21['rin'][x]],fmt='r.')
    ax2.errorbar(obstime[5]-50000,mean([100,i20e21['rinerr1'][5]]),\
                yerr=100-mean([100,i20e21['rinerr1'][5]]),\
                 uplims=True,fmt='r')
    ax2.set_xlim([4100,5400])
    ax2.set_ylim([0,100])
    ax2.set_yticklabels(('0','20','40','60','80'))
    ax2.set_ylabel('i=20,e=2.1')
    
    ax3 = plt.subplot(5,1,3)
    ax3.set_xticklabels('',visible=False)
    x = [0,3,4,5,6]
    ax3.errorbar(obstime[x]-50000,i50e30['rin'][x],yerr=[i50e30['rin'][x]-\
        i50e30['rinerr1'][x],i50e30['rinerr2'][x]-i50e30['rin'][x]],fmt='r.')
    ax3.errorbar(obstime[1:3]-50000,[mean([100,i50e30['rinerr1'][1]]),\
                mean([100,i50e30['rinerr1'][2]])],yerr=[100-\
                 mean([100,i50e30['rinerr1'][1]]),100-\
            mean([100,i50e30['rinerr1'][2]])],uplims=True,fmt='r.',visible=False)
    ax3.set_xlim([4100,5400])
    ax3.set_ylim([0,100])
    ax3.set_yticklabels(('0','20','40','60','80'))
    ax3.set_ylabel('i=50,e=3.0')
    
    ax4 = plt.subplot(5,1,4)
    ax4.set_xticklabels('',visible=False)
    x = [0,3,4,6]
    ax4.errorbar(obstime[x]-50000,i50e21['rin'][x],yerr=[i50e21['rin'][x]-\
        i50e21['rinerr1'][x],i50e21['rinerr2'][x]-i50e21['rin'][x]],fmt='r.')
    ery = array([mean([100,i50e21['rinerr1'][1]]),\
                 mean([100,i50e21['rinerr1'][2]]),\
            mean([100,i50e21['rinerr1'][5]])])
    limy = 100-ery
    ax4.errorbar(obstime[[1,2,5]]-50000,ery,yerr=limy,\
                 uplims=True,fmt='r.',visible=False)
    ax4.set_xlim([4100,5400])
    ax4.set_ylim([0,100])
    ax4.set_yticklabels(('0','20','40','60','80'))
    ax4.set_ylabel('i=50,e=2.1')
    
    ax5 = plt.subplot(5,1,5)
    ax5.set_xticklabels('',visible=False)
    for i in range(7):
        ax5.plot(np.array([obstime[i],obstime[i]])-50000,\
                 [-1000,7000],'k--',linewidth=1)
    ax5.errorbar(pca['time']-50000,pca['rate']/1000,pca['error']/1000,fmt='.')
    ax5.set_xlim([4100,5400])
    ax5.set_ylim([0,7])
    ax5.set_yticklabels(('0','1','2','3','4','5','6'))
    ax5.set_ylabel('PCA Rate (1000 Hz)')
    ax5.set_xlabel('Observation Time (MJD-50000)')
    fig.subplots_adjust(hspace=.0)

#Make dnorm vs. lc plot
#Make Rin vs Time plot
def dnormlc():
    fig = figure(1)
    plt.clf()
    ax1 = plt.subplot(5,1,1)
    ax1.set_xticklabels('',visible=False)
    plt.title('$R_{in}$ vs. Observation Time (Diskbb)')
    rin,rinerr1,rinerr2 = dnormtorg(i20e30,20)
    ax1.errorbar(obstime[[0,1,2,3,6]]-50000,rin,yerr=[rin-rinerr1,rinerr2-rin],fmt='r.')
    ax1.set_xlim([4100,5400])
##    ax1.set_ylim([0,1])
    ax1.set_ylabel('i=20,e=3.0')
    ax1.set_yticklabels(('.2','.4','.6','.8','1.0','1.2','1.4','1.6','1.8',\
                         '2.0'))
    
    ax2 = plt.subplot(5,1,2)
    ax2.set_xticklabels('',visible=False)
    rin,rinerr1,rinerr2 = dnormtorg(i20e21,20)
    ax2.errorbar(obstime[[0,1,2,3,6]]-50000,rin,yerr=[rin-rinerr1,rinerr2-rin],fmt='r.')
    ax2.set_xlim([4100,5400])
    ax2.set_ylim([.2,2.0])
    ax2.set_yticklabels(('.2','.4','.6','.8','1.0','1.2','1.4','1.6','1.8'))
    ax2.set_ylabel('i=20,e=2.1')
    
    ax3 = plt.subplot(5,1,3)
    ax3.set_xticklabels('',visible=False)
    rin,rinerr1,rinerr2 = dnormtorg(i50e30,50)
    ax3.errorbar(obstime[[0,1,2,3,4,6]]-50000,rin,yerr=[rin-rinerr1,rinerr2-rin],fmt='r.')
    ax3.set_xlim([4100,5400])
    ax3.set_yticklabels(('5','10','15','20','25','30','35'))
    ax3.set_ylabel('i=50,e=3.0')
    
    ax4 = plt.subplot(5,1,4)
    ax4.set_xticklabels('',visible=False)
    rin,rinerr1,rinerr2 = dnormtorg(i50e21,50)
    ax4.errorbar(obstime[[0,1,2,3,4,6]]-50000,rin,yerr=[rin-rinerr1,rinerr2-rin],fmt='r.')
    ax4.set_xlim([4100,5400])
    ax4.set_yticklabels(('5','10','15','20','25','30','35'))
    ax4.set_ylabel('i=50,e=2.1')
    
    ax5 = plt.subplot(5,1,5)
    ax5.set_xticklabels('',visible=False)
    for i in range(7):
        ax5.plot(np.array([obstime[i],obstime[i]])-50000,\
                 [-1000,7000],'k--',linewidth=1)
    ax5.errorbar(pca['time']-50000,pca['rate']/1000,pca['error']/1000,fmt='.')
    ax5.set_xlim([4100,5400])
    ax5.set_ylim([0,7])
    ax5.set_yticklabels(('0','1','2','3','4','5','6'))
    ax5.set_ylabel('PCA Rate (1000 Hz)')
    ax5.set_xlabel('Observation Time (MJD-50000)')
    fig.subplots_adjust(hspace=.0)
    
def ftest(d1,d2,ind=None):
    if ind is None:
        ind = range(size(d1[0]))
    dof1 = size(d1[0][ind])-10 #DoF for kdblur model is N-10 (10 free parameters)
    dof2 = size(d2[0][ind])-9 #DoF for null model is N-9 (9 free parameters)
    var1 = sum((d1[2][ind]-d1[4][ind])**2/d1[3][ind]**2)/dof1
    var2 = sum((d2[2][ind]-d2[4][ind])**2/d1[3][ind]**2)/dof2
    f = var1/var2
    pvalue = scipy.stats.f.cdf(f,dof1,dof2)
    print 'Var1: '+str(var1)
    print 'Var2: '+str(var2)
    print 'DoF1: '+str(dof1)
    print 'DoF2: '+str(dof2)
    print 'f: '+str(f)
    print 'Pvalue: '+str(pvalue)
##    if isnan(pvalue) == True:
##        pdb.set_trace()
    return pvalue, f

#Read in qdp data
def qdpread(filename):
    f = open(filename,'r')
    f = f.readlines()
    sz = size(f)
    noind = where(array(f)=='NO NO NO NO NO\n')[0]
    #Convert to float array
    for i in arange(3,noind):
        a = array(f[i].split(' '))
        #pdb.set_trace()
        f[i] = a.astype('float')
    for i in arange(noind+1,sz):
        a = array(f[i].split(' '))
        f[i] = a.astype('float')
    xrt = f[3:noind]
    pca = f[noind+1:]
    return np.transpose(xrt), np.transpose(pca)

#Script for Quentin
def qscript():
    xrt1, pca1 = qdpread('WithDiskbbi20e30.qdp') #Read more complicated model
    xrt2, pca2 = qdpread('i20e30Kdblur.qdp') #Read less complicated model
    ftest(xrt1,xrt2) #Perform f-test using only XRT data
    ftest(np.concatenate((xrt1,pca1),axis=1),concatenate((xrt2,pca2),axis=1)) #All data

#Carry out f-test on all fits
def fanalysis():
    #Initialize array
    i20e30p = np.zeros((4,7))
    i20e30f = np.zeros((4,7))
    i20e21p = np.zeros((4,7))
    i20e21f = np.zeros((4,7))
    i50e30p = np.zeros((4,7))
    i50e30f = np.zeros((4,7))
    i50e21p = np.zeros((4,7))
    i50e21f = np.zeros((4,7))

    
    os.chdir('/Users/rallured/GX339-4/')
    #Loop through 
    for s in range(7):
        #Read in i20e30 data
        xrt1, pca1 = qdpread('spectrum'+str(s+1)+'/i20e30Kdblur.qdp')
        comb1 = np.concatenate((xrt1,pca1),axis=1)
        xrt2, pca2 = qdpread('spectrum'+str(s+1)+'/NoKdblur.qdp')
        comb2 = np.concatenate((xrt2,pca2),axis=1)
        ind = where(logical_and(xrt1[0]<8.,xrt1[0]>5.))
        cind = where(logical_and(comb1[0]<8.,comb1[0]>5.))
        
        #Calculate f statistics
        i20e30p[0,s], i20e30f[0,s] = ftest(comb1,comb2) #Combined full energy range
        i20e30p[1,s], i20e30f[1,s] = ftest(comb1,comb2,cind) #Combined 5-8 keV
        i20e30p[2,s], i20e30f[2,s] = ftest(xrt1,xrt2) #XRT full energy range
        i20e30p[3,s], i20e30f[3,s] = ftest(xrt1,xrt2,ind) #XRT 5-8 keV

        #Read in i20e21 data
        xrt1, pca1 = qdpread('spectrum'+str(s+1)+'/i20e21Kdblur.qdp')
        comb1 = np.concatenate((xrt1,pca1),axis=1)
##        xrt2, pca2 = qdpread('spectrum'+str(s+1)+'/i20e21NoKdblur.qdp')
##        comb2 = np.concatenate((xrt2,pca2),axis=1)
        ind = where(logical_and(xrt1[0]<8.,xrt1[0]>5.))
        cind = where(logical_and(comb1[0]<8.,comb1[0]>5.))
        #Calculate f statistics
        i20e21p[0,s], i20e21f[0,s] = ftest(comb1,comb2) #Combined full energy range
        i20e21p[1,s], i20e21f[1,s] = ftest(comb1,comb2,cind) #Combined 5-8 keV
        i20e21p[2,s], i20e21f[2,s] = ftest(xrt1,xrt2) #XRT full energy range
        i20e21p[3,s], i20e21f[3,s] = ftest(xrt1,xrt2,ind) #XRT 5-8 keV

        #Read in i50e30 data
        xrt1, pca1 = qdpread('spectrum'+str(s+1)+'/i50e30Kdblur.qdp')
        comb1 = np.concatenate((xrt1,pca1),axis=1)
##        xrt2, pca2 = qdpread('spectrum'+str(s+1)+'/i50e30NoKdblur.qdp')
##        comb2 = np.concatenate((xrt2,pca2),axis=1)
        ind = where(logical_and(xrt1[0]<8.,xrt1[0]>5.))
        cind = where(logical_and(comb1[0]<8.,comb1[0]>5.))
        #Calculate f statistics
        i50e30p[0,s], i50e30f[0,s] = ftest(comb1,comb2) #Combined full energy range
        i50e30p[1,s], i50e30f[1,s] = ftest(comb1,comb2,cind) #Combined 5-8 keV
        i50e30p[2,s], i50e30f[2,s] = ftest(xrt1,xrt2) #XRT full energy range
        i50e30p[3,s], i50e30f[3,s] = ftest(xrt1,xrt2,ind) #XRT 5-8 keV

        #Read in i50e21 data
        xrt1, pca1 = qdpread('spectrum'+str(s+1)+'/i50e21Kdblur.qdp')
        comb1 = np.concatenate((xrt1,pca1),axis=1)
##        xrt2, pca2 = qdpread('spectrum'+str(s+1)+'/i50e21NoKdblur.qdp')
##        comb2 = np.concatenate((xrt2,pca2),axis=1)
        ind = where(logical_and(xrt1[0]<8.,xrt1[0]>5.))
        cind = where(logical_and(comb1[0]<8.,comb1[0]>5.))
        #Calculate f statistics
        i50e21p[0,s], i50e21f[0,s] = ftest(comb1,comb2) #Combined full energy range
        i50e21p[1,s], i50e21f[1,s] = ftest(comb1,comb2,cind) #Combined 5-8 keV
        i50e21p[2,s], i50e21f[2,s] = ftest(xrt1,xrt2) #XRT full energy range
        i50e21p[3,s], i50e21f[3,s] = ftest(xrt1,xrt2,ind) #XRT 5-8 keV
        
    pvalues = reshape([i20e30p,i20e21p,i50e30p,i50e21p],(4,4,7))
    fvalues = reshape([i20e30f,i20e21f,i50e30f,i50e21f],(4,4,7))

    return pvalues, fvalues

#Create simulated spectra
def simspectra():
    os.chdir('/Users/rallured/GX339-4/spectrum7/i20e30sim/')
    print 'ok'
    #Load in XRT spectrum and background
    xrtspecf = pyfits.open('combine_rebin.pha')
    xrtspec = xrtspecf[1].data
    xrtbackf = pyfits.open('combine.bak')
    xrtback = xrtbackf[1].data

    #Load in PCA spectrum and background
    pcaspecf = pyfits.open('pca_sys.pha')
    pcaspec = pcaspecf[1].data
    pcabackf = pyfits.open('pca_back.pha')
    pcaback = pcabackf[1].data
    os.chdir('/Users/rallured/GX339-4/spectrum7/i20e30sim/')

    #Initialize new spectra arrays
    nxrtspec = zeros(size(xrtspec['channel']))
    nxrtback = copy(nxrtspec)
    npcaspec = zeros(size(pcaspec['channel']))
    npcaback = copy(npcaspec)
    
    #Create 500 simulated spectra using Poisson counting statistics
    for s in range(500):
        #Create new spectra
        for i in range(size(xrtspec['channel'])):
            nxrtspec[i] = np.random.poisson(lam=xrtspec['counts'][i])
            nxrtback[i] = np.random.poisson(lam=xrtback['counts'][i])
        for i in range(size(pcaspec['channel'])):
            npcaspec[i] = np.random.poisson(lam=pcaspec['counts'][i])
            npcaback[i] = np.random.poisson(lam=pcaback['counts'][i])
        
        #Create new tables
        xrtspecf[1].data.field('counts')[:] = nxrtspec
        xrtspecf.writeto('fake_combine_rebin'+str(s+1)+'.pha')
        
        xrtbackf[1].data.field('counts')[:] = nxrtback
        xrtbackf.writeto('fake_combine_back'+str(s+1)+'.pha')

        pcaspecf[1].data.field('counts')[:] = npcaspec
        pcaspecf.writeto('fake_pca_sys'+str(s+1)+'.pha')

        pcabackf[1].data.field('counts')[:] = npcaback
        pcabackf.writeto('fake_pca_back'+str(s+1)+'.pha')

#Parse show files and create XSPEC setup models
def parseshow(inp,out):
    print 'ok-1'
    fin = open(inp,'r')
    print 'ok-2'
    s = fin.readlines()
    print 'ok-3'
    mod = open(inp.split('/')[0]+'/model02_reflionx.xcm','r')
    mod = mod.readlines()
    print mod
    fout = open(out,'w')
    print 'ok0'
    fout.writelines(mod[0:9]) #Set up model header
    print 'ok1'
    
    if s[7][12:16] == 'disk':#Normal model (diskbb)
        for p in range(12):
            fout.write(s[6+p][43:54]+mod[9+p][20:])
            print 'ok' + str(p+1)
        fout.write('=9\n')
        for p in arange(13,16):
            fout.write(s[6+p][43:54]+mod[9+p][20:])
        fout.write(s[23][43:54]+mod[25][20:])
        fout.writelines(mod[26:])
    else:#Diskbb omitted
        for p in range(10):
            fout.write(s[6+p][43:54]+mod[9+p][20:])
        fout.write('=7\n')
        for p in arange(11,14):
            fout.write(s[6+p][43:54]+mod[9+p][20:])
        fout.write(s[21][43:54]+mod[23][20:])
        fout.writelines(mod[24:])

#Parse show files and create XSPEC setup models
def parsenull(inp,out):
    print 'ok-1'
    fin = open(inp,'r')
    print 'ok-2'
    s = fin.readlines()
    print 'ok-3'
    mod = open(inp.split('/')[0]+'/model02_reflionx_null.xcm','r')
    mod = mod.readlines()
    print mod
    fout = open(out,'w')
    print 'ok0'
    fout.writelines(mod[0:9]) #Set up model header
    print 'ok1'
    
    if s[7][12:16] == 'disk':#Normal model (diskbb)
        for p in range(8):
            fout.write(s[6+p][43:54]+mod[9+p][20:])
            print 'ok' + str(p+1)
        fout.write('=5\n')
        for p in arange(9,12):
            fout.write(s[6+p][43:54]+mod[9+p][20:])
        fout.write(s[19][43:54]+mod[21][20:])
        fout.writelines(mod[22:])
    else:#Diskbb omitted
        for p in range(6):
            fout.write(s[6+p][43:54]+mod[9+p][20:])
        fout.write('=3\n')
        for p in arange(7,10):
            fout.write(s[6+p][43:54]+mod[9+p][20:])
        fout.write(s[17][43:54]+mod[19][20:])
        fout.writelines(mod[20:])
    fout.close()
    print fout
    print os.getcwd()

#Loop through spectra and create all XSPEC models
def parseallspectra():
    os.chdir('/Users/rallured/GX339-4/')
    files = glob.glob('spectrum*/NullFit.txt')

    for f in files:
##        out = '/Users/rallured/GX339-4/'+f.split('/')[0]+'/i'+\
##              f.split('wi')[1][:5]+'Model2.xcm'
        out = '/Users/rallured/GX339-4/'+f.split('/')[0]+'/NullFit.xcm'
        parsenull(f,out)

#Make p value plots for f-test
def pvalueplots():
    pvalues, fvalues = fanalysis()
    plt.clf()
    plt.hist((pvalues[:,0,0],pvalues[:,0,1],pvalues[:,0,2],pvalues[:,0,3],\
          pvalues[:,0,4],pvalues[:,0,5],pvalues[:,0,6]),histtype='bar',\
         label=['Spectrum1','Spectrum2','Spectrum3','Spectrum4','Spectrum5',\
                'Spectrum6','Spectrum7'])
    plt.title('F-test P-Values (XRT+PCA)')
    plt.xlabel('P-Value')
    plt.ylabel('Occurrences')
    plt.legend(loc='upper left')
    plt.savefig('PValBoth.eps')

    plt.clf()
    plt.hist((pvalues[:,2,0],pvalues[:,2,1],pvalues[:,2,2],pvalues[:,2,3],\
          pvalues[:,2,4],pvalues[:,2,5],pvalues[:,2,6]),histtype='bar',\
         label=['Spectrum1','Spectrum2','Spectrum3','Spectrum4','Spectrum5',\
                'Spectrum6','Spectrum7'])
    plt.title('F-test P-Values (XRT Only)')
    plt.xlabel('P-Value')
    plt.ylabel('Occurrences')
    plt.legend(loc='upper left')
    plt.savefig('PValXRT.eps')

    plt.clf()
    pvalues = pvalues[:,0,:]
    

#Diskbb normalization to radius in units of Rg
def dnormtorg(d,i):
    #Correction factors
    cor = 1.7**2 * .412
    if i==20:
        ind = [0,1,2,3,4,5]
    elif i==50:
        ind = [0,1,2,3,4,5]
    else:
        ind = [0,1,2,3,4,5]
        i = 50
    i = i*pi/180.
    rin = cor * sqrt(d['dnorm'][ind]/cos(i)) * (4./5.)
    rinerr1 = cor * sqrt(d['dnormerr1'][ind]/cos(i)) * (4./5.)
    rinerr2 = cor * sqrt(d['dnormerr2'][ind]/cos(i)) * (4./5.)
    g_c2 = 7.42E-31
    solar = 1.99E30
    rg = g_c2*solar*10#5.8/(sin(i)**3)
    return rin/rg, rinerr1/rg, rinerr2/rg

#Make plot of dnorm vs. pca flux
def dnormvsflux():
    timelag = 0.
    x = np.zeros(6)
    c = 0
    for i in [0,1,2,3,4,6]:
        diff = pca['time']-(obstime[i]-timelag)
        close = where(abs(diff)==min(abs(diff)))[0]
        if pca['time'][close] > obstime[i]-timelag:
            close = close-1
        x[c] = close
        c = c+1
    x = x.astype('int')
    rate = pca['rate']
    rate = rate[x]
    rerr = pca['error']
    rerr = rerr[x]
    print x
            
##        slope = (pca['rate'][close+1]-pca['rate'][close])\
##                /(pca['time'][close+1]-pca['time'][close])
##        lum[i] = pca['rate'][close] + slope*\
##                 (obstime[i]-timelag-pca['time'][close])
##        diff = pca['time']-(obstime[i])
##        close = where(abs(diff)==min(abs(diff)))[0]
##        if pca['time'][close] > obstime[i]:
##            close = close-1
##        slope = (pca['rate'][close+1]-pca['rate'][close])\
##                /(pca['time'][close+1]-pca['time'][close])
##        lum2[i] = pca['rate'][close] + slope*\
##                 (obstime[i]-pca['time'][close]
    
    plt.clf()
##    rin, rinerr1, rinerr2 = dnormtorg(i50e30,50)
##    plt.errorbar(rate,rin,xerr=rerr\
##       ,yerr=[rin-rinerr1,rinerr2-rin],fmt='b.',ecolor='b')
    rin, rinerr1, rinerr2 = dnormtorg(i50e21,50)
    plt.errorbar(rate,rin,xerr=rerr\
        ,yerr=[rin-rinerr1,rinerr2-rin],fmt='r.',ecolor='r')
##    noblur = np.genfromtxt('/Users/rallured/GX339-4/NoBlur.txt',\
##        names='nh,nherr1,nherr2,tin,tinerr1,tinerr2,dnorm,dnormerr1,\
##dnormerr2,phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
##,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
##    rin, rinerr1, rinerr2 = dnormtorg(noblur,0)
##    plt.errorbar(rate,rin,xerr=rerr\
##        ,yerr=[rin-rinerr1,rinerr2-rin],fmt='r.',ecolor='r')    
    plt.xlabel('PCA Rate (cts/sec)')
    plt.ylabel('$R_{in}$ ($GM/c^2$)')
    plt.title('$R_{in}$ vs. PCA Rate (Diskbb)')

noblur = np.genfromtxt('/Users/rallured/GX339-4/NoBlur.txt',\
        names='nh,nherr1,nherr2,tin,tinerr1,tinerr2,dnorm,dnormerr1,\
        dnormerr2,phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
        ,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')

f = array([2.0035e-9,7.0981e-10,1.8365e-9,1.9048e-9,2.8518e-9,8.847e-9])
ftol = 4*pi*(2.469E22)**2
f = f*ftol
lum = f/1.3e39

def lumplot(d,k=True,i=20,e=3.0):
    f = array([2.0035e-9,7.0981e-10,1.8365e-9,1.9048e-9,2.8518e-9,8.847e-9])
    ftol = 4*pi*(2.469E22)**2
    f = f*ftol
    ledd = 1.3E39 #eg/s
    print f/ledd
    lum = f/ledd

    
    if k==False:
        rin, rinerr1, rinerr2 = dnormtorg(d,0)
    else:
        rin, rinerr1, rinerr2 = [d['rin'],d['rinerr1'],d['rinerr2']]
        lum = append(lum,.0223)
        if i==20 and e==2.1:
            rin = append(rin,9.3381)
            rinerr1 = append(rinerr1,1.235)
            rinerr2 = append(rinerr2,24.9847)
        elif i==20 and e==3.0:
            rin = append(rin,16.284)
            rinerr1 = append(rinerr1,11.111)
            rinerr2 = append(rinerr2,27.9365)
        elif i==20 and e==4.0:
            rin = append(rin,19.8428)
            rinerr1 = append(rinerr1,15.7114)
            rinerr2 = append(rinerr2,32.656)
        elif i==50 and e==2.1:
            rin = append(rin,32.1942)
            rinerr1 = append(rinerr1,8.63477)
            rinerr2 = append(rinerr2,63.9984)
        elif i==50 and e==3.0:
            rin = append(rin,39.961)
            rinerr1 = append(rinerr1,24.9975)
            rinerr2 = append(rinerr2,63.9911)
        elif i==50 and e==4.0:
            rin = append(rin,44.3858)
            rinerr1 = append(rinerr1,32.5946)
            rinerr2 = append(rinerr2,64.0101)
    
##    f = array([1.912e-9,6.796e-10,1.752e-9,1.834e-9,2.775e-9,8.501e-9]) #erg/s/cm^2
##    f = array([3.319e-9,1.7e-9,3.324e-9,2.835e-9,3.548e-9,1.149e-8]) #erg/s/cm^2
    matplotlib.rcParams['legend.fontsize']=14
##    plt.clf()
    plt.semilogy(lum,rin,visible=False)
    plt.errorbar(lum,rin,yerr=[(rin-rinerr1),(rinerr2-rin)],\
                 fmt='k.',ecolor='r',\
                 label='Our measurements',markersize=3)
    plt.title('i='+str(i)+',e='+str(e))
    plt.xlabel('Source Luminosity ($L_{Edd}$)')
    plt.ylabel('Inner Radius ($GM/c^2$)')
##    plt.errorbar(array([.14,.14,.46,1.33,2.0,3.25,3.25])/100,\
##                 [75,190,2.9,3.6,13.3,4.0,6.33],\
##                 yerr=[[10,90,.7,1,6,.5,.08],[10,710,2.1,1.4,6.4,.5,1.04]]\
##                 ,lolims=[True,False,False,False,False,False,\
##                False],fmt='b.',ecolor='r')
    plt.errorbar([.14/100],[75],yerr=[10],uplims=True,fmt='k.',ecolor='g',\
                 label='Tomsick et al. (2009)',markersize=3)
    plt.errorbar([.14/100,.02],[190,13.3],yerr=[[90,6.],[710.,6.4]],fmt='k.',ecolor='y',\
                 label='Shidatsu et al. (2011)',markersize=3)
    plt.errorbar([.0046,.0133],[2.9,3.6],yerr=[[.7,1.],[2.1,1.4]],fmt='k.',\
                 ecolor='b',label='Tomsick et al. (2008)',markersize=3)
    plt.errorbar([.0325],[4],yerr=[.5],fmt='k.',ecolor='orange',label='Miller et al. (2006)')
    plt.errorbar([.0325,.0325],[6.33,24],yerr=[[.08,3],[1.04,1]],fmt='k.',\
                 ecolor='c',label='Done & Diaz Trigo (2010)',markersize=3)
    plt.errorbar([.0325],[2.8],yerr=[.1],fmt='k.',markersize=3,ecolor='m',\
                 label='Reis et al. (2008)')
##    plt.legend(loc='upper right')
    plt.text(.034,22.,'i=60 (fixed)',fontsize=14)
    plt.text(.034,6.4,'i=22 (fit)',fontsize=14)
    #Calculate Rin constant chisq
    guess = arange(1.236,200,.01)
    cfit = fitconst(rin,(rin-rinerr1)/1.645,(rinerr2-rin)/1.645,guess)
    guess = guess[where(isnan(cfit)==False)]
    cfit = cfit[where(isnan(cfit)==False)]
##    plt.clf()
##    plt.plot(guess,cfit)
    print min(cfit)*5
    print guess[where(cfit == min(cfit))]
    pvalue = 1-scipy.stats.chi2.cdf(min(cfit)*5,5)
    print pvalue
    matplotlib.rcParams['legend.fontsize']=16

#Make radius plots for all fitting rounds
def radplots(k=True):
    clf()
    subplot(3,2,1)
    lumplot(i20e21,k=k,i=20,e=2.1)
    subplot(3,2,3)
    lumplot(i20e30,k=k,i=20,e=3.0)
    subplot(3,2,5)
    lumplot(i20e40,k=k,i=20,e=4.0)
    subplot(3,2,2)
    lumplot(i50e21,k=k,i=50,e=2.1)
    subplot(3,2,4te)
    lumplot(i50e30,k=k,i=50,e=3.0)
    subplot(3,2,6)
    lumplot(i50e40,k=k,i=50,e=4.0)

#Is Fe constant?
def feconst(noblur):
    guess = arange(.8,2.0,.001)
##    fe,feerr1,feerr2 = append(noblur['fe'],.966736),append(noblur['feerr1'],.809944),append(noblur['feerr2'],1.34188)
##    fe,feerr1,feerr2 = append(noblur['fe'],.952622),append(noblur['feerr1'],.785177),append(noblur['feerr2'],1.27517)
##    fe,feerr1,feerr2 = append(noblur['fe'],.935724),append(noblur['feerr1'],.781455),append(noblur['feerr2'],1.25789)
##    fe,feerr1,feerr2 = append(noblur['fe'],1.14008),append(noblur['feerr1'],.897471),append(noblur['feerr2'],1.40861)
##    fe,feerr1,feerr2 = append(noblur['fe'],1.1416),append(noblur['feerr1'],.893475),append(noblur['feerr2'],1.35043)
    fe,feerr1,feerr2 = append(noblur['fe'],1.16925),append(noblur['feerr1'],.89276),append(noblur['feerr2'],1.45106)
    cfit = fitconst(fe,(fe-feerr1)/1.645,\
                    (feerr2-fe)/1.645,guess)
    print min(cfit)*5
    print guess[where(cfit == min(cfit))]
    print 1-scipy.stats.chi2.cdf(min(cfit)*5,5)

#Make dbp plot
def dbp():
    dbp = np.genfromtxt('/Users/rallured/GX339-4/spectrum4/DBP.txt')
    plt.clf()
    plt.plot(arange(1,8),dbp[0,:],'h',label='i=20,e=3.0',markersize=10)
    plt.plot(arange(1,8),dbp[1,:],'s',label='i=20,e=2.1',markersize=10)
    plt.plot(arange(1,8),dbp[2,:],'p',label='i=50,e=3.0',markersize=10)
    plt.plot(arange(1,8),dbp[3,:],'*',label='i=50,e=2.1',markersize=10)
    plt.legend(loc='upper left')
    plt.xlim([0,8])
    plt.ylabel('P-Value')
    plt.xlabel('Spectrum $\#$')
    plt.title('Diskbb F-test P-Values (XRT+PCA)')
    plt.savefig('DiskPVal.eps')

    pvalues, fvalues = fanalysis()
    pvalues = pvalues[:,0,:]
    plt.clf()
    plt.plot(arange(1,8),pvalues[0,:],'h',label='i=20,e=3.0',markersize=10)
    plt.plot(arange(1,8),pvalues[1,:],'s',label='i=20,e=2.1',markersize=10)
    plt.plot(arange(1,8),pvalues[2,:],'p',label='i=50,e=3.0',markersize=10)
    plt.plot(arange(1,8),pvalues[3,:],'*',label='i=50,e=2.1',markersize=10)
    plt.legend(loc='upper left')
    plt.xlim([0,8])
    plt.ylabel('P-Value')
    plt.xlabel('Spectrum $\#$')
    plt.title('Kdblur F-test P-Values (XRT+PCA)')
    plt.savefig('BlurPVal.eps')
    
    return dbp

#Make plot using inner radii from nH fixed, i free fits
def incldet():
    noblur = np.genfromtxt('/Users/rallured/GX339-4/NoBlur.txt',\
        names='nh,nherr1,nherr2,tin,tinerr1,tinerr2,dnorm,dnormerr1,\
        dnormerr2,phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
        ,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')
    rin, rinerr1, rinerr2 = dnormtorg(noblur,0)
    os.chdir('/Users/rallured/GX339-4/Results120606/')
    d = np.zeros((6,30))
    i = 0
    for f in glob.glob('*/*Param.txt'):
        t = np.genfromtxt(f).flatten()
        if size(t)==30:
            d[i,:] = t
            i = i+1
    d = np.transpose(d)

    plt.clf()
    plt.loglog(d[9],rin,visible=False)
    plt.errorbar(d[9],rin,xerr=[d[9]-d[10],d[11]-d[9]],yerr=[rin-rinerr1,\
            rinerr2-rin],fmt='b.',ecolor='r')

#Make nh dnorm cor plot
def nhdnorm():
    plt.clf()
    plt.loglog(noblur['dnorm'],noblur['nh'],'.',visible=False)
    plt.errorbar(noblur['dnorm'],noblur['nh'],xerr=[noblur['dnorm']\
        -noblur['dnormerr1'],noblur['dnormerr2']-noblur['dnorm']],\
        yerr=[noblur['nh']-noblur['nherr1'],noblur['nherr2']-noblur['nh']],\
                 fmt='.')

nhfixed = np.genfromtxt('/Users/rallured/GX339-4/nHFixed.txt',\
        names='tin,tinerr1,tinerr2,dnorm,dnormerr1,\
        dnormerr2,phoind,phoinderr1,phoinderr2,pnorm,pnormerr1\
        ,pnormerr2,fe,feerr1,feerr2,xi,xierr1,xierr2,rnorm,rnormerr1,rnormerr2')

#Calculate average column densities for all 6 geometries
def avgnH():
    nh = []
    guess = arange(.6,1.2,.001)
    cfit = fitconst(i20e21['nh'],(i20e21['nh']-i20e21['nherr1'])/1.645,\
                    (i20e21['nherr2']-i20e21['nh'])/1.645,guess)
    nh.append(guess[where(cfit == min(cfit))])
    cfit = fitconst(i20e30['nh'],(i20e30['nh']-i20e30['nherr1'])/1.645,\
                    (i20e30['nherr2']-i20e30['nh'])/1.645,guess)
    nh.append(guess[where(cfit == min(cfit))])
    cfit = fitconst(i20e40['nh'],(i20e40['nh']-i20e40['nherr1'])/1.645,\
                    (i20e40['nherr2']-i20e40['nh'])/1.645,guess)
    nh.append(guess[where(cfit == min(cfit))])
    cfit = fitconst(i50e21['nh'],(i50e21['nh']-i50e21['nherr1'])/1.645,\
                    (i50e21['nherr2']-i50e21['nh'])/1.645,guess)
    nh.append(guess[where(cfit == min(cfit))])
    cfit = fitconst(i50e30['nh'],(i50e30['nh']-i50e30['nherr1'])/1.645,\
                    (i50e30['nherr2']-i50e30['nh'])/1.645,guess)
    nh.append(guess[where(cfit == min(cfit))])
    cfit = fitconst(i50e40['nh'],(i50e40['nh']-i50e40['nherr1'])/1.645,\
                    (i50e40['nherr2']-i50e40['nh'])/1.645,guess)
    nh.append(guess[where(cfit == min(cfit))])
    return nh

#Compute a mean correlation coefficient assuming normal statistics
def weightedrho(xvec,yvec,yerr1,yerr2):
    #Convert to 3 sigma
    yerr1 = (yvec-yerr1)*3/1.645
    yerr2 = (yerr2-yvec)*3/1.645
    scanlength = 5 #Number of data points to scan in each direction

    #Create y loop vectors for each observation
    yscan = zeros((6,scanlength*2+1))
    for i in range(6):
        ylow = yvec[i] - arange(scanlength,-1,-1)*yerr1[i]/scanlength
        yhigh = yvec[i] + arange(1,scanlength+1)*yerr2[i]/scanlength
        yscan[i] = concatenate((ylow,yhigh))

    #Create weight vector (identical for each observation)
    weights = zeros(scanlength*2+1)
    p = arange(-3.,3.+3./scanlength,3./scanlength)
    for i in range(scanlength*2+1):
        weights[i] = scipy.stats.norm.cdf(p[i]+1.5/scanlength)-\
                     scipy.stats.norm.cdf(p[i]-1.5/scanlength)

    #Perform integral
    integral = 0.
    integral2 = 0.
    norm = 0.
    for x0 in range(2*scanlength+1):
        for x1 in range(2*scanlength+1):
            for x2 in range(2*scanlength+1):
                for x3 in range(2*scanlength+1):
                    for x4 in range(2*scanlength+1):
                        for x5 in range(2*scanlength+1):
                            vec = array([yscan[0,x0],yscan[1,x1],yscan[2,x2],\
                            yscan[3,x3],yscan[4,x4],yscan[5,x5]])
                            coeff = scipy.stats.spearmanr(lum,vec)[0]
                            integral = integral + weights[x0]*weights[x1]*\
    weights[x2]*weights[x3]*weights[x4]*weights[x5]*coeff
                            integral2 = integral2 + weights[x0]*weights[x1]*\
    weights[x2]*weights[x3]*weights[x4]*weights[x5]*coeff**2
                            norm = norm + weights[x0]*weights[x1]*\
    weights[x2]*weights[x3]*weights[x4]*weights[x5]
                            print str(x0)+'\t'+str(x1)+'\t'+str(x2)+'\t'+str(x3)\
                                  +'\t'+str(x4)+'\t'+str(x5)
                            sys.stdout.flush()

    pdb.set_trace()

    return [integral,integral2,norm]

#Generate normal value on one side of distribution
def getnorm(mean,sigma,pos=True):
    x = scipy.stats.norm.rvs(loc=mean,scale=sigma)
    if x >= 0 and pos==True:
        return x
    elif x < 0 and pos==False:
        return x
    else:
        return getnorm(mean,sigma,pos=pos)

#Simulate coefficient distribution
def coeffdist(xvec,yvec,yerr1,yerr2):
    #Convert to 1 sigma
    yerr1 = (yvec-yerr1)/1.645
    yerr2 = (yerr2-yvec)/1.645
    yerr = array([yerr1,yerr2])

    length = 10000
    yscan = zeros((6,length))
    #Generate rin values, first generate postive/negative vector, then loop
    #through and call getnorm to generate distributions
    for i in range(6):
        posneg = scipy.stats.binom.rvs(1,.5,size=length)
        #Handle pegged at lower limit
        if yerr1[i] == 0:
            posneg = repeat(1,length)
        for j in range(length):
            yscan[i,j] = yvec[i] + getnorm(0,yerr[posneg[j],i],pos=posneg[j])

    coeff = zeros(length)
    #Compute coefficients
    for j in range(length):
        coeff[j] = scipy.stats.pearsonr(lum,yscan[:,j])[0]

    return coeff

#Create correlation plots
def corplots(cor=None):
    #i20 kdblur plot
    if cor is None:
        wi20e21r = coeffdist(append(lum,.0223),append(i20e21['rin'],9.3381),\
                append(i20e21['rinerr1'],1.235),append(i20e21['rinerr2'],24.9847))
        wi20e30r = coeffdist(append(lum,.0223),append(i20e30['rin'],16.284),\
                append(i20e30['rinerr1'],11.111),append(i20e30['rinerr2'],27.9365))
        wi20e40r = coeffdist(append(lum,.0223),append(i20e40['rin'],19.8428),\
                append(i20e40['rinerr1'],15.7114),append(i20e40['rinerr2'],32.656))
##        wi20e21d = coeffdist(lum,i20e21['dnorm'],i20e21['dnormerr1'],i20e21['dnormerr2'])
        wi20e21d = coeffdist(lum,*dnormtorg(i20e21,20))
##        wi20e30d = coeffdist(lum,i20e30['dnorm'],i20e30['dnormerr1'],i20e30['dnormerr2'])
        wi20e30d = coeffdist(lum,*dnormtorg(i20e30,20))
##        wi20e40d = coeffdist(lum,i20e40['dnorm'],i20e40['dnormerr1'],i20e40['dnormerr2'])
        wi20e40d = coeffdist(lum,*dnormtorg(i20e40,20))
        wi50e21r = coeffdist(append(lum,.0223),append(i50e21['rin'],32.1942),\
                append(i50e21['rinerr1'],8.63477),append(i50e21['rinerr2'],63.9984))
        wi50e30r = coeffdist(append(lum,.0223),append(i50e30['rin'],39.961),\
                append(i50e30['rinerr1'],24.9975),append(i50e30['rinerr2'],63.9911))
        wi50e40r = coeffdist(append(lum,.0223),append(i50e40['rin'],44.3858),\
                append(i50e40['rinerr1'],32.5946),append(i50e40['rinerr2'],64.0101))
        wi50e21d = coeffdist(lum,*dnormtorg(i50e21,50))
        wi50e30d = coeffdist(lum,*dnormtorg(i50e30,50))
        wi50e40d = coeffdist(lum,*dnormtorg(i50e40,50))
        fi20e21 = coeffdist(lum,i20e21['fe'],i20e21['feerr1'],i20e21['feerr2'])
        fi20e30 = coeffdist(lum,i20e30['fe'],i20e30['feerr1'],i20e30['feerr2'])
        fi20e40 = coeffdist(lum,i20e40['fe'],i20e40['feerr1'],i20e40['feerr2'])
        fi50e21 = coeffdist(lum,i50e21['fe'],i50e21['feerr1'],i50e21['feerr2'])
        fi50e30 = coeffdist(lum,i50e30['fe'],i50e30['feerr1'],i50e30['feerr2'])
        fi50e40 = coeffdist(lum,i50e40['fe'],i50e40['feerr1'],i50e40['feerr2'])
        
    else:
        wi20e21r = cor[0]
        wi20e30r = cor[1]
        wi20e40r = cor[2]
        wi20e21d = cor[3]
        wi20e30d = cor[4]
        wi20e40d = cor[5]
        wi50e21r = cor[6]
        wi50e30r = cor[7]
        wi50e40r = cor[8]
        wi50e21d = cor[9]
        wi50e30d = cor[10]
        wi50e40d = cor[11]
        fi20e21 = cor[12]
        fi20e30 = cor[13]
        fi20e40 = cor[14]
        fi50e21 = cor[15]
        fi50e30 = cor[16]
        fi50e40 = cor[17]

    #Print mean correlations
    if cor is not None:
        for i in range(18):
            print str(mean(cor[i])) + ' ' + str(std(cor[i]))

    clf()
    subplot(3,2,1)
    hist(wi20e21r,bins=100,histtype='step',label='e=2.1')
    hist(wi20e30r,bins=100,histtype='step',label='e=3.0')
    hist(wi20e40r,bins=100,histtype='step',label='e=4.0')
    legend(loc='upper right')
    title('Kdblur Correlation Distribution (i=20)\n\n')
    xlabel('Pearson Correlation Coefficient')
    ylabel('# of Occurrences')
##    savefig('/Users/rallured/GX339-4/130409Results/i20KdblurCor.eps')

    #i20 dnorm plot
    subplot(3,2,3)
    hist(wi20e21d,bins=100,histtype='step',label='e=2.1')
    hist(wi20e30d,bins=100,histtype='step',label='e=3.0')
    hist(wi20e40d,bins=100,histtype='step',label='e=4.0')
##    legend(loc='upper left')
    title('Diskbb Correlation Distribution (i=20)\n\n')
    xlabel('Pearson Correlation Coefficient')
    ylabel('# of Occurrences')
##    savefig('/Users/rallured/GX339-4/130409Results/i20DnormCor.eps')

    #i50 kdblur plot
    subplot(3,2,2)
    hist(wi50e21r,bins=100,histtype='step',label='e=2.1')
    hist(wi50e30r,bins=100,histtype='step',label='e=3.0')
    hist(wi50e40r,bins=100,histtype='step',label='e=4.0')
##    legend(loc='upper left')
    title('Kdblur Correlation Distribution (i=50)\n\n')
    xlabel('Pearson Correlation Coefficient')
    ylabel('# of Occurrences')
##    savefig('/Users/rallured/GX339-4/130409Results/i50KdblurCor.eps')

    #i50 dnorm plot
    subplot(3,2,4)
    hist(wi50e21d,bins=100,histtype='step',label='e=2.1')
    hist(wi50e30d,bins=100,histtype='step',label='e=3.0')
    hist(wi50e40d,bins=100,histtype='step',label='e=4.0')
##    legend(loc='upper left')
    title('Diskbb Correlation Distribution (i=50)\n\n')
    xlabel('Pearson Correlation Coefficient')
    ylabel('# of Occurrences')
##    savefig('/Users/rallured/GX339-4/130409Results/RinCorPlots.eps')

    #i20 fe plot
    subplot(3,2,5)
    hist(fi20e21,bins=100,histtype='step',label='e=2.1')
    hist(fi20e30,bins=100,histtype='step',label='e=3.0')
    hist(fi20e40,bins=100,histtype='step',label='e=4.0')
##    legend(loc='upper left')
    title('Fe Correlation Distribution (i=20)\n\n')
    xlabel('Pearson Correlation Coefficient')
    ylabel('# of Occurrences')
##    savefig('/Users/rallured/GX339-4/130409Results/RinCorPlots.eps')

    #i50 fe plot
    subplot(3,2,6)
    hist(fi50e21,bins=100,histtype='step',label='e=2.1')
    hist(fi50e30,bins=100,histtype='step',label='e=3.0')
    hist(fi50e40,bins=100,histtype='step',label='e=4.0')
##    legend(loc='upper left')
    title('Fe Correlation Distribution (i=50)\n\n')
    xlabel('Pearson Correlation Coefficient')
    ylabel('# of Occurrences')
    savefig('/Users/rallured/GX339-4/130409Results/CorPlots.eps')

    #Return correlation calculations
    return (wi20e21r,wi20e30r,wi20e40r,wi20e21d,wi20e30d,wi20e40d,wi50e21r,\
            wi50e30r,wi50e40r,wi50e21d,wi50e30d,wi50e40d,fi20e21,fi20e30,fi20e40,\
            fi50e21,fi50e30,fi50e40)
