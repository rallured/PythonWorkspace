from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
import glob,pdb
from plotting import *
from zernikemod import locateimage

#Function to make 3d contour plot of file
def makeplot(filename,zmin,zmax,fig):
    d = genfromtxt(filename)
    d = -d * 1000.
    cx,cy,rad = locateimage(d)
    x = arange(128.)
    x,y = meshgrid(x,x)
    r = sqrt((x-cx)**2+(y-cy)**2)
    d[r>rad*.97] = nan
    d = d - nanmin(d)
    fig.clf()
    ax = fig.add_subplot(111,projection='3d')
    mycontour(d,nlev=250,fmt='%i')
    ax.set_zlim((zmin,zmax))
    fig.show()
