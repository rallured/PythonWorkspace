import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

from numpy import *
import time

mod = SourceModule("""
__global__ void doublify(double *a)
{
    int idx = threadIdx.x + threadIdx.y*4;
    a[idx] *= 2;
}
""")

testmod = SourceModule("""
__global__ void testsumff(float *a_in,float *a_out)
{
    //Declare shared memory
    extern __shared__ float sdataf[];

    //Load global memory into shared memory
    int tid = threadIdx.x;
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    sdataf[tid] = a_in[i];// * a_in[i]/2 + a_in[i]/10;
    __syncthreads();

    //Performed reversed, threadId-based indexed reduction
    for (int s=blockDim.x/2; s>0; s = s/2)
    {
        if (tid < s)
        {
            sdataf[tid] = sdataf[tid] + sdataf[tid+s];
        }
        __syncthreads();
    }

    //0th block element goes to global output array
    if (tid==0) a_out[blockIdx.x] = sdataf[0];
}

__global__ void testsum(double *a_in,double *a_out)
{
    //Declare shared memory
    extern __shared__ double sdata[];

    //Load global memory into shared memory
    int tid = threadIdx.x;
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    sdata[tid] = a_in[i];// * a_in[i]/2 + a_in[i]/10;
    __syncthreads();

    //Performed reversed, threadId-based indexed reduction
    for (int s=blockDim.x/2; s>0; s = s/2)
    {
        if (tid < s)
        {
            sdata[tid] = sdata[tid] + sdata[tid+s];
        }
        __syncthreads();
    }

    //0th block element goes to global output array
    if (tid==0) a_out[blockIdx.x] = sdata[0];
}

__global__ void flops(float *a_in,float b)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    a_in[i] = a_in[i]*a_in[i]/b + a_in[i]/10;
}

__global__ void dflops(double *a_in,double b)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    a_in[i] = a_in[i]*a_in[i]/b + a_in[i]/10;
}

__global__ void testsumd(double *a_in,double *a_out)
{
    //Declare shared memory
    extern __shared__ double sdata[];

    //Load global memory into shared memory
    int tid = threadIdx.x;
    int i = blockIdx.x*blockDim.x*2 + threadIdx.x;
    sdata[tid] = a_in[i] + a_in[i+blockDim.x];// * a_in[i]/2 + a_in[i]/10;
    __syncthreads();

    //Performed reversed, threadId-based indexed reduction
    for (int s=blockDim.x/2; s>32; s = s/2)
    {
        if (tid < s)
        {
            sdata[tid] = sdata[tid] + sdata[tid+s];
        }
        __syncthreads();
    }

    if (tid < 32)
    {
        sdata[tid] += sdata[tid+32];
        sdata[tid] += sdata[tid+16];
        sdata[tid] += sdata[tid+8];
        sdata[tid] += sdata[tid+4];
        sdata[tid] += sdata[tid+2];
        sdata[tid] += sdata[tid+1];
    }

    //0th block element goes to global output array
    if (tid==0) a_out[blockIdx.x] = sdata[0];
}

__global__ void testsumf(float *a_in,float *a_out)
{
    //Declare shared memory
    extern __shared__ float sdataf[];

    //Load global memory into shared memory
    int tid = threadIdx.x;
    int i = blockIdx.x*blockDim.x*2 + threadIdx.x;
    sdataf[tid] = a_in[i] + a_in[i+blockDim.x];// * a_in[i]/2 + a_in[i]/10;
    __syncthreads();

    //Performed reversed, threadId-based indexed reduction
    for (int s=blockDim.x/2; s>=32; s = s/2)
    {
        if (tid < s)
        {
            sdataf[tid] = sdataf[tid] + sdataf[tid+s];
        }
        __syncthreads();
    }

    if (tid < 32)
    {
        sdataf[tid] += sdataf[tid+32];
        sdataf[tid] += sdataf[tid+16];
        sdataf[tid] += sdataf[tid+8];
        sdataf[tid] += sdataf[tid+4];
        sdataf[tid] += sdataf[tid+2];
        sdataf[tid] += sdataf[tid+1];
    }

    //0th block element goes to global output array
    if (tid==0) a_out[blockIdx.x] = sdataf[0];
}
""")

##a = random.randn(4,4)
##
##a_gpu = cuda.mem_alloc(a.nbytes)
##
##cuda.memcpy_htod(a_gpu,a)
##
##func = mod.get_function("doublify")
##func(a_gpu, block=(4,4,1))
##
##a_doubled = empty_like(a)
##cuda.memcpy_dtoh(a_doubled,a_gpu)
##print a_doubled
##print a

#Do an intensive computation both in serial and then in parallel
#Serial
arr = random.randn(1024*1024*10)
arr = arr.astype(float32)
tstart = time.time()
print sum(arr*arr/2+arr/10)
print time.time() - tstart

#Parallel
#Launch Kernel to reduce each block to find the minimum
t0 = cuda.Event()
t1 = cuda.Event()
t2 = cuda.Event()
tstart = time.time()
arr_out = zeros(size(arr)/1024)
arr_out = arr_out.astype(float32)
d_in = cuda.mem_alloc(arr.nbytes)
d_out = cuda.mem_alloc(arr_out.nbytes)
cuda.memcpy_htod(d_in,arr)
cuda.memcpy_htod(d_out,arr_out)
testsum = testmod.get_function("testsumf")
dflops = testmod.get_function("flops")
tstart = time.time()
t0.record()
dflops(d_in,float32(2),block=(1024,1,1),grid=(size(arr_out),1,1))
testsum(d_in,d_out,block=(1024,1,1),grid=(size(arr_out)/2,1,1),shared=1024*4)
testsum(d_out,d_out,block=(1024,1,1),grid=(size(arr_out)/1024/2,1,1),shared=1024*4)
t1.record()
t1.synchronize()
print 'Kernel execution: ' + str(time.time()-tstart)
cuda.memcpy_dtoh(arr_out,d_out)
print 'Memcpy: ' + str(time.time()-tstart)
t2.record()
t2.synchronize()
print 'GPU execution: ' + str(t1.time_since(t0))
print 'GPU Memcpy: ' + str(t2.time_since(t1))

print sum(arr_out[:10])
print time.time()-tstart
