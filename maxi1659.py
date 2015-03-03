from os import *
import numpy as np
import scipy as sp
from string import *
import glob
import matplotlib.pyplot as plt
import matplotlib

def getdata():
    data = np.genfromtxt('dutchhardratio.out')
    return data

def loopthrough(hard,soft,obsid):
    for path in glob.glob('95*'):
        if path != '95358-01-01-00':
            chdir(path)
            data = getdata()
            print data
            soft.append(data[0])
            hard.append(data[1])
            obsid.append(path)
            chdir('..')
    return 0

def hardratios():
    # Initialize data vectors
    hard = []
    soft = []
    obsid = []
    # Go to parent directory
    chdir('/Users/ryanallured/RXTE/MAXI_J1659-152/results/')
    # Loop through P95108
    chdir('P95108')
    loopthrough(hard,soft,obsid)
    # Loop through P95118
    chdir('../P95118')
    loopthrough(hard,soft,obsid)
    # Loop through P95358
    chdir('../P95358')
    loopthrough(hard,soft,obsid)
    # Write data to file
    chdir('/Users/ryanallured/RXTE/MAXI_J1659-152/results/')
    np.savetxt('DutchRatios.txt',np.transpose([hard,soft,obsid]),fmt='%s   %s   %s')
    return 0

def readratios():
    # Read soft and hard fluxes from file
    chdir('/Users/ryanallured/RXTE/MAXI_J1659-152/results/')
    data = np.genfromtxt('DutchRatios.txt',names='soft,hard,obsid',dtype=None)
    return data

def plotratios(flux):
    flux
    plt.interactive(True)
    plt.hold(False)
    plt.hist(flux['hard']/flux['soft'])

def getstddata():
    p95108 = '/Users/ryanallured/RXTE/MAXI_J1659-152/results/P95108/'
    p95118 = '/Users/ryanallured/RXTE/MAXI_J1659-152/results/P95118/'
    p95358 = '/Users/ryanallured/RXTE/MAXI_J1659-152/results/P95358/'

    chdir(p95108)
    obs = glob.glob('95108-??-??-??')
    for path in obs:
        chdir(p95108+path)
        if 'data' not in locals():
            data = [np.genfromtxt('stdmodel.out',dtype=None)]
            names = [path]
        elif 'data' in locals():
            data.append(np.genfromtxt('stdmodel.out',dtype=None))
            names.append(path)

    chdir(p95118)
    obs = glob.glob('95118-??-??-??')
    for path in obs:
        chdir(p95118+path)
        data.append(np.genfromtxt('stdmodel.out',dtype=None))
        names.append(path)

    chdir(p95358)
    obs = glob.glob('95358-??-??-??')
    for path in obs:
        if path != '95358-01-01-00':
            chdir(p95358+path)
            data.append(np.genfromtxt('stdmodel.out',dtype=None))
            names.append(path)
            
    return [data,names]

def getallstddata():
    chdir('/Users/rallured/MAXI_J1659-152/results/StdResults/')
    spectral = np.genfromtxt("SpectralResults.txt",names='mflux,mfluxerr1,mfluxerr2,chisq,nh,tin, \
        tinerr1,tinerr2,phoind,phoinderr1,phoinderr2,pnorm,pnormerr1,pnormerr2,dnorm, \
        dnormerr1,dnormerr2,unabsflux,diskflux,sflux,hflux')
    timing = np.genfromtxt('TimingResults.txt',names='rms,qpofreq,qpopow',skip_header=1)
    names = np.genfromtxt('Names.txt',dtype=str)
    tstart = np.genfromtxt('Tstart.txt')

    return [spectral,timing,names,tstart]

def classifystates(spectral,timing):
    disk = spectral['diskflux']/spectral['unabsflux']

    #Thermal states
    thermal = np.logical_and(disk > .75,timing['qpopow'] < .005)
    thermal = np.logical_and(thermal,timing['rms'] < .075)

    #Hard states
    hard = np.logical_and(disk < .2,spectral['phoind'] > 1.4)
    hard = np.logical_and(hard,spectral['phoind'] < 2.1)
    hard = np.logical_and(hard,timing['rms'] > .1)

    #SPL states
    spl = np.logical_and(spectral['phoind'] > 2.4,np.logical_and(timing['rms'] < .15,\
        np.logical_or(np.logical_and(disk < .8,timing['qpopow'] > .01),disk < .5)))

    #Int states
    inter = np.logical_not(np.logical_or(thermal,np.logical_or(spl,hard)))

    return [thermal,hard,spl,inter]

def makeplots(spectral,timing,states,tstart):
    #Define secondary spectral parameters
    disk = spectral['diskflux']/spectral['unabsflux']
    hardr = spectral['hflux']/spectral['sflux']
    inters = states[3]
    inters[0] = False
    inters[1] = False
    inters[3] = False

    matplotlib.rcParams['axes.titlesize']=18
    matplotlib.rcParams['axes.labelsize']=18
    matplotlib.rcParams['xtick.labelsize']=18
    matplotlib.rcParams['ytick.labelsize']=18

    #Change to plot directory
    chdir('/Users/rallured/MAXI_J1659-152/results/StdResults/Plots/')

    #Make HID
    plt.clf()
    plt.hold(True)
    plt.plot(hardr[states[0]],spectral['mflux'][states[0]],'r.',label='Thermal',markersize=20)
    plt.plot(hardr[states[1]],spectral['mflux'][states[1]],'b.',label='Hard',markersize=20)
    plt.plot(hardr[states[2]],spectral['mflux'][states[2]],'g.',label='SPL',markersize=20)
    plt.plot(hardr[states[3]],spectral['mflux'][states[3]],'y.',label='Intermediate',markersize=20)
    plt.title('Intensity vs. Hardness')
    plt.xlabel('Hardness Ratio (6.3keV-10.5keV)/(3.8keV-7.3keV)')
    plt.ylabel(r'Model Flux (cts/s/cm$^2$)')
    plt.savefig('HID.eps')

    #Make lightcurve
    plt.clf()
    merr1 = spectral['mflux']-spectral['mfluxerr1']
    merr2 = spectral['mfluxerr2']-spectral['mflux']
    plt.plot(tstart[states[0]],spectral['mflux'][states[0]],'r.',label='Thermal',markersize=20)
    plt.plot(tstart[states[1]],spectral['mflux'][states[1]],'b.',label='Hard',markersize=20)
    plt.plot(tstart[states[2]],spectral['mflux'][states[2]],'g.',label='SPL',markersize=20)
    plt.plot(tstart[states[3]],spectral['mflux'][states[3]],'y.',label='Intermediate',markersize=20)
    plt.legend()
    plt.title('Lightcurve for MAXI J1659-152')
    plt.xlabel('Days Since Sep. 25 2010')
    plt.ylabel(r'Model Flux (cts/s/cm$^2$')
    plt.savefig('Lightcurve.eps')

    #Make power-law flux vs. disk flux
    plt.clf()
    powerlaw = spectral['unabsflux'] - spectral['diskflux']
    plt.plot(spectral['diskflux'][states[0]],powerlaw[states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(spectral['diskflux'][states[1]],powerlaw[states[1]],'.b',label='Hard',markersize=20)
    plt.plot(spectral['diskflux'][states[2]],powerlaw[states[2]],'.g',label='SPL',markersize=20)
    plt.plot(spectral['diskflux'][inters],powerlaw[inters],'.y',label='Intermediate',markersize=20)
    plt.title('Power Law Flux vs. Disk Flux')
    plt.xlabel(r'Disk Flux (cts/s/cm$^2$)')
    plt.ylabel(r'Power Law Flux (cts/s/cm$^2$)')
    plt.savefig('PlawvsDflux.eps')

    #Make disk fraction plot
    plt.clf()
    diskok = np.logical_and(spectral['dnorm'] < 1000,spectral['dnorm'] > 1)
    thermal = np.logical_and(diskok,states[0])
    hard = np.logical_and(diskok,states[1])
    spl = np.logical_and(diskok,states[2])
    intermediate = np.logical_and(diskok,states[3])
    plt.plot(tstart[thermal],disk[thermal],'.r',label='Thermal',markersize=20)
    plt.plot(tstart[hard],disk[hard],'.b',label='Hard',markersize=20)
    plt.plot(tstart[spl],disk[spl],'.g',label='SPL',markersize=20)
    plt.plot(tstart[intermediate],disk[intermediate],'.y',label='Intermediate',markersize=20)
    plt.title('Disk Fraction vs. Observation Time')
    plt.xlabel('Days Since Sep. 25 2010')
    plt.ylabel('Disk Fraction')
    plt.legend(loc='lower center')
    plt.savefig('Disk.eps')

    plt.clf()
    #Make Photon Index vs. Hardness
    plt.subplot(3,1,1)
    plt.plot(hardr[states[0]],spectral['phoind'][states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(hardr[states[1]],spectral['phoind'][states[1]],'.b',label='Hard',markersize=20)
    plt.plot(hardr[states[2]],spectral['phoind'][states[2]],'.g',label='SPL',markersize=20)
    plt.plot(hardr[states[3]],spectral['phoind'][states[3]],'.y',label='Intermediate',markersize=20)
    plt.yticks(np.arange(1.7,2.6,.2))
    plt.xticks(np.arange(.2,.65,.05),())
    plt.title('Spectral Parameters vs. Hardness Ratio')
    plt.ylabel(r'$\Gamma$')
    #Make Disk fraction vs. hardness
    plt.subplot(3,1,2)
    plt.plot(hardr[states[0]],disk[states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(hardr[states[1]],disk[states[1]],'.b',label='Hard',markersize=20)
    plt.plot(hardr[states[2]],disk[states[2]],'.g',label='SPL',markersize=20)
    plt.plot(hardr[inters],disk[inters],'.y',label='Intermediate',markersize=20)
    plt.yticks(np.arange(.0,.7,.2))
    plt.ylabel('f')
    plt.xticks(np.arange(.2,.65,.05),())
    #Make RMS power vs hardness
    plt.subplot(3,1,3)
    plt.plot(hardr[states[0]],timing['rms'][states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(hardr[states[1]],timing['rms'][states[1]],'.b',label='Hard',markersize=20)
    plt.plot(hardr[states[2]],timing['rms'][states[2]],'.g',label='SPL',markersize=20)
    plt.plot(hardr[states[3]],timing['rms'][states[3]],'.y',label='Intermediate',markersize=20)
    plt.yticks(np.arange(0.,.25,.05))
    plt.subplots_adjust(hspace=0.)
    plt.ylabel('r')
    plt.xlabel('Hardness Ratio (6.3keV-10.5keV)/(3.8keV-7.3keV)')
    plt.savefig('SpectralvsHard.eps')

    #Make Photon Index vs. Tstart
    plt.clf()
    plt.plot(tstart[states[0]],spectral['phoind'][states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(tstart[states[1]],spectral['phoind'][states[1]],'.b',label='Hard',markersize=20)
    plt.plot(tstart[states[2]],spectral['phoind'][states[2]],'.g',label='SPL',markersize=20)
    plt.plot(tstart[states[3]],spectral['phoind'][states[3]],'.y',label='Intermediate',markersize=20)
    plt.title('Photon Index vs. Time')
    plt.ylabel('Photon Index')
    plt.xlabel('Days Since Sep. 25 2010')
    plt.savefig('PhoIndvsTime.eps')

    #Make RMS power vs Tstart
    plt.clf()
    plt.plot(tstart[states[0]],timing['rms'][states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(tstart[states[1]],timing['rms'][states[1]],'.b',label='Hard',markersize=20)
    plt.plot(tstart[states[2]],timing['rms'][states[2]],'.g',label='SPL',markersize=20)
    plt.plot(tstart[states[3]],timing['rms'][states[3]],'.y',label='Intermediate',markersize=20)
    plt.title('RMS Power vs. Time')
    plt.ylabel('RMS Power')
    plt.xlabel('Days Since Sep. 25 2010')
    plt.savefig('RMSvsTime.eps')

    #Make QPO power vs Tstart
    plt.clf()
    plt.plot(tstart[states[0]],timing['qpopow'][states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(tstart[states[1]],timing['qpopow'][states[1]],'.b',label='Hard',markersize=20)
    plt.plot(tstart[states[2]],timing['qpopow'][states[2]],'.g',label='SPL',markersize=20)
    plt.plot(tstart[states[3]],timing['qpopow'][states[3]],'.y',label='Intermediate',markersize=20)
    plt.title('QPO Power vs. Time')
    plt.ylabel('QPO Power')
    plt.xlabel('Days Since Sep. 25 2010')
##    plt.subplots_adjust(left=.15)
    plt.savefig('QPOvsTime.eps')
##    plt.subplots_adjust(left=.125)

    #Make D_norm vs tstart
    plt.clf()
    yerr = [np.sqrt(spectral['dnorm']-spectral['dnormerr1']),\
            np.sqrt(spectral['dnormerr2']-spectral['dnorm'])]
    plt.semilogy(tstart,np.sqrt(spectral['dnorm']),visible=False)
    plt.errorbar(tstart,np.sqrt(spectral['dnorm']),fmt='.',yerr=yerr)
    plt.title(r'$\sqrt{Disk\/Normalization}$')
    plt.xlabel('Days Since Sep. 25 2010')
    plt.savefig('DnormvsTime.eps')

    matplotlib.rcParams['axes.titlesize']=18
    matplotlib.rcParams['axes.labelsize']=16
    matplotlib.rcParams['xtick.labelsize']=14
    matplotlib.rcParams['ytick.labelsize']=14
