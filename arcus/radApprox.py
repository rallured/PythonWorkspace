#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from numba import cuda
import numba
from accelerate.cuda.blas import Blas
import math
import time,pdb
from accelerate.cuda import rand
prng = rand.PRNG(rndtype=rand.PRNG.XORWOW)

MAX32 = np.uint32(0xffffffff)

@cuda.jit("(uint64[::1], uint64)", device=True)
def cuda_xorshift(states, id):
    x = states[id]
    x ^= x >> 12
    x ^= x << 25
    x ^= x >> 27
    states[id] = x
    return numba.uint64(x) * numba.uint64(2685821657736338717)


@cuda.jit("float32(uint64[::1], uint64)", device=True)
def cuda_xorshift_float(states, id):
    return numba.float32(numba.float32(MAX32 & cuda_xorshift(states, id)) \
                      / numba.float32(MAX32))



def createGratingArray(Ng=1e5,snapto=7):
    """
    Create an array of grooves that are discretized to
    a given snap-to length.
    Snap-to is achieved through the use of the round
    function and is defined as an order of magnitude.
    snapto=7 would be a 0.1 nm snapto
    
    """
    r = np.linspace(12e3-5.,12e3+5,1e5)
    dt = 160e-6/12e3 #fan angle groove to groove
    theta = np.arange(Ng)*dt
    theta = theta - np.mean(theta)

    for i in range(Ng):
        #Compute groove coordinates
        x = r*np.sin(theta[i])
        z = r*np.cos(theta[i])-12e3
        #Apply snapto
        x = round(x,snapto)
        z = round(z,snapto)
        y = np.zeros(len(r))
        #Rotate about x axis at center of grating
        y = z*np.sin(1.5*np.pi/180)
        z = z*np.cos(1.5*np.pi/180)
        #Determine distance vector to focus
        d2 = np.sqrt((12e3+z)**2+x**2+y**2)
        #Determine distance vector to observation point(s)
        
#Trying CUDA reduction
@cuda.jit('void(complex128[:],double,double,double,'
          'double,double,double,double)')
def sections(outArr,x0,y0,z0,period,zmin,zmax,offset):
    """
    Divide the grating up into discrete sections with
    straight grooves.
    In order for each thread to identify with a specific
    groove, same number of grooves assumed for each section.
    Launch a Kernel for each section with option to add
    small constant offset to groove positions.
    """
    #Which groove?
    thisgroove = cuda.grid(1)
    #How many grooves?
    groovenum = cuda.blockDim.x*1024

    t = cuda.threadIdx.x + cuda.blockIdx.x*cuda.blockDim.x

    #Load portion of array into memory?
##    tempArr = inpArr[groovesec*1024:(groovesec+1)*1024]
    #Create random array?
    

    #Perform sum
    out = np.double(0.)
    x = numba.float64((thisgroove - groovenum/2) * period)
    for i in range(int((zmax-zmin)/.01)):
        #Create x,y,z position
##        fanangle = 160e-9/12*(thisgroove-groovenum/2)
        z = numba.float64(zmin + i*.01)-12e3
        #Rotate to optical axis reference frame
        y = z * math.sin(1.5*np.pi/180)
        z = z * math.cos(1.5*np.pi/180)

        #Determine distance vector to focus
        d1 = math.sqrt((12e3+z)**2+x**2+y**2)

        #Determine distance vector to observation point
        d2 = math.sqrt((12e3+z-z0)**2+(x-x0)**2+(y-y0)**2)

        #Determine phase factor
        dl = (d1-d2)*2*np.pi/1e-6
        out += math.cos(dl) + 1j*math.sin(dl)

    #Save to output index groovenum + groovesection
    outArr[thisgroove] = out/d2

#Compute intensity at a field position
def intSections(x0,y0,z0,gnum=10*1024,secnum=5):
    #Set up device array
    out = np.zeros(gnum)*1j
    outg = cuda.to_device(out)

    #call it!
    stream = cuda.stream()
    tstart = time.time()
    zsec = np.linspace(12e3-25.,12e3+25.,secnum+1)
    outt = 0.*1j
    offsets = np.random.uniform(low=0.,high=160e-9,size=secnum)
    for i in range(secnum):
        d = (np.mean(zsec[i:i+2]))/12e3*160e-6
        sections[gnum/1024,1024](outg,x0,y0,z0,d,zsec[i],zsec[i+1],\
                                 0.)
        stream.synchronize()
        outg.copy_to_host(out)
        outt += out.sum()
##    pdb.set_trace()
    return abs(outt.sum())**2

#Trying CUDA reduction
@cuda.jit('void(complex128[:],double,double,double,double)')
def redEx(outArr,x0,y0,z0,hubdist):
    """
    Each thread loads an array, maybe does some operations
    on it, and then sums the array.
    """
    #Round precision
    roundPrec = 10.**7
    #Which groove?
    thisgroove = cuda.grid(1)
    #How many grooves?
    groovenum = cuda.blockDim.x*1024

    t = cuda.threadIdx.x + cuda.blockIdx.x*cuda.blockDim.x

    #Load portion of array into memory?
##    tempArr = inpArr[groovesec*1024:(groovesec+1)*1024]


    #Perform sum
    out = np.double(0.)
    fanangle = 160e-6/hubdist*(thisgroove-groovenum/2)
    for i in range(1024*5):
        #Create x,y,z position
        zt = np.double((i-2560)*.01+hubdist)
        #zt = z[i]
        x = zt * math.sin(fanangle)
        zt = zt * math.cos(fanangle) - hubdist
        #Rounding operation goes here
##        z = math.floor(z*roundPrec+.5)/roundPrec
##        x = math.floor(x*roundPrec+.5)/roundPrec
        #Rotate to optical axis reference frame
        y = zt * math.sin(1.5*np.pi/180)
        zt = zt * math.cos(1.5*np.pi/180)

        #Determine distance vector to focus
        d1 = math.sqrt((12e3+zt)**2+x**2+y**2)

        #Determine distance vector to observation point
        d2 = math.sqrt((12e3+zt-z0)**2+(x-x0)**2+(y-y0)**2)

        #Determine phase factor
        dl = (d1-d2)*2*np.pi/1e-6
        out += math.cos(dl) + 1j*math.sin(dl)

    #Save to output index groovenum + groovesection
    outArr[thisgroove] = out/d2

#Need to write Kernel to reduce the groove array
#This will cut down on memory transfer time, regardless
#of whether reduction Kernel is fully optimized
@cuda.jit('void(complex128[:],complex128[:])')
def reduceGrooves(inpArr,outArr):
    #Load the array into shared memory
    sarr = cuda.shared.array(shape=0,dtype=numba.complex128)
    i = cuda.grid(1)
    tid = cuda.threadIdx.x
    if i < len(inpArr):
        sarr[tid] = inpArr[i]
    else:
        sarr[tid] = 0.*1j
    cuda.syncthreads()

    #Do reduction in shared memory
    s = cuda.blockDim.x/2
    while s>0:
        if tid < s:
            sarr[tid] += sarr[tid+s]
        cuda.syncthreads()
        s = s/2

    #Write output to global memory
    if tid==0:
        outArr[cuda.blockIdx.x] = sarr[0]
    
    return None

#Compute intensity at a field position
def intensity(x0,y0,z0,gnum=10*1024,hubdist=12e3):
    #Set up device array
    out = np.zeros(gnum)*1j
    outg = cuda.to_device(out)

##    #Random z
##    zr = np.random.uniform(low=12e3-25.,high=12e3+25.,size=1024*5)

    #call it!
    stream = cuda.stream()
    tstart = time.time()
    redEx[gnum/1024,1024](outg,x0,y0,z0,hubdist)
    stream.synchronize()
    print time.time()-tstart

    #Get output
    tstart = time.time()
    outg.copy_to_host(out)
    print time.time()-tstart

    return abs(out.sum())**2

#Evaluate operation of Kernel with NVVP
x = intensity(0,0,0)
