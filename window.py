from numpy import *
from matplotlib.pyplot import *
import os
import glob
import pdb
import sys

os.chdir('/Users/ryanallured/Documents/Research/GEMS/Detector/Window/ThermalTest')
data = genfromtxt('ThermalData.txt')
timevec = data[0]
tempvec = data[1]
presvec = data[2]
tempsig = data[3]
del data
os.chdir('../LeakTest')

def analyzeRGA(file):
    data = transpose(genfromtxt(file,dtype='string',skip_header=2))
    pressure = zeros(45)
    fdata = zeros(450)
    for i in range(450):
        fdata[i] = float(data[1][i].split('"')[1])
    print fdata
    for i in range(45):
        pressure[i] = mean(fdata[i*10:(i+1)*10])
    return pressure

def analyzeTrend(file):
    data = transpose(genfromtxt(file,dtype='string',skip_header=2))
    fdata = zeros(size(data[1]))
    for i in range(size(data[1])):
        fdata[i] = float(data[1][i].split('"')[1])
    return fdata

##os.chdir('/Users/ryanallured/Documents/Research/GEMS/Detector/Window/ThermalTest/')
##
##d1 = transpose(genfromtxt('ThermalTestResults1.txt'))
##d2 = transpose(genfromtxt('ThermalTest_2012_3_10_20_53.txt'))
##d2[0] = d2[0] + d1[0][-1]+80*60
##
##timevec = append(d1[0],d2[0])
##tempvec = append(d1[1],d2[1])
##presvec = append(d1[2],d2[2])
##
##
##os.chdir('/Users/ryanallured/Documents/Research/GEMS/Detector/Window/ThermalTest/Splice')
##
##for file in glob.glob('Thermal*'):
##    temp = transpose(genfromtxt(file))
##    temp[0] = temp[0] + 95*60**2 + 37*60 + 66679.74
##    if sum(isnan(temp[0])) != 0:
##           pdb.set_trace()
##    timevec = append(timevec,temp[0])
##    tempvec = append(tempvec,temp[1])
##    presvec = append(presvec,temp[2])
##    
##os.chdir('/Users/ryanallured/Documents/Research/GEMS/Detector/Window/ThermalTest/Splice2')
##
##for file in glob.glob('Thermal*'):
##    temp = transpose(genfromtxt(file))
##    temp[0] = temp[0] + 265*60**2 + 13*60 + 66679.74
##    if sum(isnan(temp[0])) != 0:
##           pdb.set_trace()
##    timevec = append(timevec,temp[0])
##    tempvec = append(tempvec,temp[1])
##    presvec = append(presvec,temp[2])
##
###Sort
##timevec, ind = unique(timevec,return_index=True)
##tempvec = tempvec[ind]
##presvec = presvec[ind]
##
##timevec2 = zeros(21346)
##tempvec2 = zeros(21346)
##presvec2 = zeros(21346)
##tempsig = zeros(21346)
###Rebin data into minutes
##for m in range(21346):
##    ind = where(logical_and(timevec < (m+1)*60,timevec > m*60))
##    if size(ind) != 0:
##        timevec2[m] = mean(timevec[ind])
##        presvec2[m] = mean(presvec[ind])
##        tempvec2[m] = mean(tempvec[ind])
##        tempsig[m] = std(tempvec[ind])
##        print m
##        sys.stdout.flush()
##del timevec
##del tempvec
##del presvec
