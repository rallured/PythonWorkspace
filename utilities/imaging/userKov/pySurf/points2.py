import numpy as np
from matplotlib import pyplot as plt
import math
from scipy import interpolate as ip
from numpy import ma
import os
import pdb
from copy import deepcopy
from pyGeneralRoutines.span import span

"""
Module containing functions acting on a point cloud. With the intention
of creating a class. Points are in format (Npoints,Ndim).
The class was created in pointCloud.py. In this version I am trying to maintain
the old python routines embedding them in the class. 
"""


method='linear'


"""2014/03/17 moved here all routines from EA2019_align and from OP1S07_align."""

"""2014/03/06 from v1. Start from nanovea txt files use pure python to avoid coordinates mess. Assume data are saved with y fast axis and +,+ scanning
directions. It is not guaranteed to work with different scanning directions."""

"""angpos copiata qui, che non so come si fa ad importare in questo casino di directories. 
CHANGES: angpos e' stata corretta funzionava con dati trasposti"""
    
def _angpos2(xy):
    '''given an image and some notable points, return the angular positions of points with respect to barycenter.
    The angle returned is in the range [-pi:pi]'''

    #N xy points are expressed as array N x 2
    b=(np.sum(xy,0))/(xy.shape[0])    #calculate barycenter
    deltaxy=xy-(b.repeat(xy.size/2).reshape(2,-1)).T
    r=np.sqrt(np.sum(deltaxy**2,1))
    theta=np.arctan2(deltaxy[:,1],deltaxy[:,0])
    return theta,r,b

def rotate_points(points,theta,center=(0,0)):
    """returns rotated coordinates of 2D point(s) x ([Npoints x 2]) about a center with anticlockwise angle theta in rad. If 3D points are passed, z coordinate is maintained."""
    tx,ty=center
    if (points.shape[-1]==3):
        return np.hstack([rotate_points(points[:,0:2],theta,center),points[:,2:]])
    else:
        if(points.shape[-1]!=2):
            raise Exception
    x,y=np.array(points[:,0]),np.array(points[:,1])
    cost=math.cos(theta)
    sint=math.sin(theta)
    x1=x*cost-y*sint + ty*sint - tx*(cost-1)
    y1=x*sint+y*cost - tx*sint - ty *(cost-1)
    return np.vstack((x1,y1)).T

def translate_points(x,offset=(0,0)):
    """returns translated coordinates of 2D point(s) x ([Npoints x 2]) by an offset.
    It works also on 3-D points, in that case the z is returned unchanged."""
    points=np.matrix(x) 
    #2014/08/11 
    #- return a 3d array if a 3d array was passed as argument.
    #- convert the result to array before return.
    #   The result was initially of type np.matrix to simplify operations on extracted columns,
    #   but it didn't really work. The problem is that the result is a matrix and if the
    #   user is unaware of that, all following slicing (column extractions) will be messed up
    #   (vectors will be column vectors as opposite to row vectors returned by np arrays.
    
    x,y=np.array(points[:,0]),np.array(points[:,1])
    x1=x+offset[0]
    y1=y+offset[1]
    translated=np.hstack((x1,y1))
    if points.shape[1]==3:
        translated=np.hstack([x1,y1,points[:,2]])
    return np.array(translated)

def get_points(filename,x=None,y=None,xrange=None,yrange=None,matrix=False,scale=(1.,1.,1.),center=None,skip_header=None,delimiter=','):
    """return a set of xyz points (N,3) from nanovea saved txt (matrix=False)
    or gwyddion saved matrix (matrix=True, xrange, yrange must be defined).
    Scale is used to scale the 
    Center is the position of the center of the image (before any scaling or rotation) in absolute coordinates.
        Set to (0,0) to center the coordinate system to the data."""
    #2014/04/29 added x and y as preferred arguments to xrange and yrange (to be removed).
    if skip_header==None: 
        skip=0 
    else: 
        skip=skip_header
        
    if (matrix):
        if ((xrange==None or yrange==None) and (x==None or y==None)): raise Exception
        mdata=np.genfromtxt(filename,skip_header=skip)
        nx,ny= mdata.shape
        if x==None and y==None:
            x=np.linspace(*xrange,num=nx)
            y=np.linspace(*yrange,num=ny)
        #plt.imshow(mdata, extent=(mxrange.min(), myrange.max(), mxrange.max(), myrange.min()))
        xpoints,ypoints=[xy.flatten() for xy in np.array(np.meshgrid(x,y))]
        zpoints=mdata.flatten()
        points=np.vstack([xpoints,ypoints,zpoints]).T
    else:
        points= np.genfromtxt(filename,delimiter=delimiter,skip_header=skip)
    
    if center != None:
        offset=np.hstack([((np.nanmax(points,axis=0)-np.nanmin(points,axis=0))/2)[0:2],0])-np.array([center[0],center[1],0])
        points=points-offset  #center on 0  
    if scale != None:
        points=points*scale
    return points

def save_points(points,filename,xgrid=None,ygrid=None,shape=None,matrix=False,fill_value=np.nan,**kwargs):
    """save points on a file. If matrix is true write in matrix form (in this case you have to 
    provide the values for axis). Otherwise write as points in columns."""
    #2014/08/08 default fill_value modified to nan.
    #20140420 moved filename to second argument for consistency.
    if matrix:
        if shape==None:
            assert len(xgrid.shape)==len(ygrid.shape)==1
        else:
            assert xgrid==ygrid==None
            x,y,z=np.hsplit(points,3)
            xgrid=np.linspace(x.min(),x.max(),shape[0])
            ygrid=np.linspace(y.min(),y.max(),shape[1])
        assert xgrid!=None
        assert ygrid!=None
        grid=np.vstack([g.flatten() for g in  np.meshgrid(xgrid,ygrid)]).T
        points=ip.griddata(points[:,0:2],points[:,2],grid,fill_value=fill_value,method=method)
        points=points.reshape(ygrid.size,xgrid.size)
        #points=np.hstack([grid,points[:,np.newaxis]])
    #if not, they are already in the correct format
    np.savetxt(filename,points,**kwargs)

def resample_points(tpoints,positions):
    """resample tpoints [Npoints x 3] on the points defined in positions [Mpoints x 2], or [Mpoints x 3]
    (in this case 3rd column is ignored).
    Return a [Nx x Ny , 3] vector of points. To get a (plottable) matrix of data use:
    plt.imshow(rpoints[:,2].reshape(xgrid.size,ygrid.size))."""
    assert tpoints.shape[1]==3
    z=ip.griddata(tpoints[:,0:2],tpoints[:,-1],positions[:,0:2],method=method)
    rpoints=np.hstack([positions[:,0:2],z[:,np.newaxis]])
    return rpoints

def resample_grid(tpoints,xgrid,ygrid):
    from matplotlib.mlab import griddata
    """resample tpoints [Npoints x 3] on the grid defined by two vectors xgrid [Nx] and ygrid [Ny].
    Return a [Nx * Ny , 3] vector of points, sorted in standard python order
    (x changes faster). To get a (plottable) matrix of data use:
    plt.imshow(rpoints[:,2].reshape(ygrid.size,xgrid.size))."""
    assert tpoints.shape[1]==3
    x,y=np.meshgrid(xgrid,ygrid) #this always return
    z=ip.griddata(tpoints[:,0:2],tpoints[:,2],(x,y),method=method) 
    rpoints=np.vstack([x.flatten(),y.flatten(),z.flatten()]).T   
    #2014/11/24 rpoints=np.vstack([x.T.flatten(),y.T.flatten(),z.T.flatten()]).T
    return rpoints
    
def plot_points(points,xgrid=None,ygrid=None,shape=None,resample=True,scatter=False,**kwargs):
    """resample xyz points [Nx3] to a grid whose axis xgrid and ygrid are given
    and plot it. If resample is set to False x and y positions are considered only for range, 
    but they are not used to position the z values (it works if x and y are on an exact unrotated grid,
    resampling is slower, but exact)."""
    
    #plot
    #plt.clf()
    x,y,z=np.hsplit(points,3)
    cmap=kwargs.pop('cmap','jet')
    aspect=kwargs.pop('aspect','equal')
    #pdb.set_trace()
    if scatter: #do scatterplot
        beamsize=20 #I may use this to represent beamsize (round symbol) or lateral resolution (square)
        plt.gca().set_aspect(aspect)
        plt.scatter(x, y, c=z, s=beamsize, cmap=cmap, edgecolors='None', **kwargs)
    else:
        #if not scatter plot, grid must be provided as shape or as xgrid and ygrid axis.
        if shape==None:
            if not len(xgrid.shape)==len(ygrid.shape)==1:
                return z    #skip plot
        else:
            assert xgrid==ygrid==None
            xgrid=np.linspace(x.min(),x.max(),shape[0])
            ygrid=np.linspace(y.min(),y.max(),shape[1])
        if resample:
            print "resampling..."
            z=resample_grid(points,xgrid,ygrid)[:,2]
        nx,ny=[xx.size for xx in (xgrid,ygrid)]
        xr,yr=[span(xx) for xx in (xgrid,ygrid)]
        xxg,sx=np.linspace(xr[0],xr[1],nx,retstep=1)
        yyg,sy=np.linspace(yr[0],yr[1],ny,retstep=1)
        #ranges for plot (intervals centered on xxg, yyg)
        xr=xr+np.array((-sx/2,sx/2))
        yr=yr+np.array((-sy/2,sy/2))
        
        z=z.reshape(ny,nx)
        plt.imshow(z,extent=[xr[0],xr[1],yr[0],yr[1]],interpolation='none',aspect=aspect,
            origin='lower')
    plt.xlabel('X (mm)')
    plt.ylabel('Y (mm)')
    plt.colorbar()
    plt.show()
    return z

def points_rms(points,xcrop=None,ycrop=None,mask=None,clip=None):
    """return the rms of the selected points."""
    
    """Note that if instead of deepcopy is used:
    pp=points[:]
    the clipping operation on p2 affects the values of points.
    """
    pp=deepcopy(points)
    z=pp[:,2]
    zm=np.sqrt(np.nanmean(z**2))
    if clip != None:
        if clip[0]!=None:
            z[np.where(z<=clip[0])]=zm
        if clip[1]!=None:
            z[np.where(z>=clip[1])]=zm
    return np.sqrt(np.nanmean(z**2))
    

def crop_points(points, xrange=None,yrange=None):
    """crop a xyz points [Nx3], keeping only points inside xrange and yrange defined as (min,max)."""
    if xrange != None: 
        points=points[np.where(xrange[0]<points[:,0]),:][0]
        points=points[np.where(points[:,0]<xrange[1]),:][0]
    if yrange != None: 
        points=points[np.where(yrange[0]<points[:,1]),:][0]
        points=points[np.where(points[:,1]<yrange[1]),:][0]
    return points
    
def subtract_points(p1,p2,xySecond=False):
    """Subtract second set of points after interpolation on first set coordinates.
    If xySecond is set to True data are interpolate on xy of p2 and then subtracted."""
    if not(xySecond):
        p2=resample_points(p2,p1)
    else:
        p1=resample_points(p1,p2)
    p1[:,2]=p1[:,2]-p2[:,2]
    return p1
    
    
def smooth_points(points,xywidth):
    """"""
    def gauss_kern(size, sizey=None):
        """ Returns a normalized 2D gauss kernel array for convolutions """
        size = int(size)
        if not sizey:
            sizey = size
        else:
            sizey = int(sizey)
        x, y = mgrid[-size:size+1, -sizey:sizey+1]
        g = exp(-(x**2/float(size)+y**2/float(sizey)))
        return g / g.sum()

    def blur_image(im, n, ny=None) :
        """ blurs the image by convolving with a gaussian kernel of typical
            size n. The optional keyword argument ny allows for a different
            size in the y direction.
        """
        g = gauss_kern(n, sizey=ny)
        improc = signal.convolve(im,g, mode='valid')
        return(improc)
    
    
    x,y,z=np.hsplit(points,3)
    if shape==None:
        assert len(xgrid.shape)==len(ygrid.shape)==1
    else:
        assert xgrid==ygrid==None
        xgrid=np.linspace(x.min(),x.max(),shape[0])
        ygrid=np.linspace(y.min(),y.max(),shape[1])
    if resample:
        print "resampling..."
        z=resample_grid(points,xgrid,ygrid)[:,2]
    z=z.reshape(xgrid.size,ygrid.size)
    z=blur_image(z,xwidth,ywidth)
    plt.imshow(z,extent=[ygrid.min(),ygrid.max(),xgrid.max(),xgrid.min()],aspect=aspect,interpolation='none')
    plt.xlabel('X (mm)')
    plt.ylabel('Y (mm)')
    plt.colorbar()
    plt.show()
    return z

    
class Points(np.ndarray):
    """Represent a set of points in a N-dimensional point cloud, defined as array in format (Npoints,Ndim)."""
    
    def __new__(cls, filename, *args,**kwargs):
        points=load(filename,*args,**kwargs)
        obj = np.asarray(points).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.info = getattr(obj, 'info', None)

    def __array_wrap__(self, out_arr, context=None):
        return np.ndarray.__array_wrap__(self, out_arr, context)   

