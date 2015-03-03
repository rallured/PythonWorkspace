from numpy import *
import socket
from scipy.optimize import curve_fit
import time

class picoMotors:
    def __init__(self,ipaddr):
        self.s = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        self.s.connect((ipaddr,23))
##        self.s.send('SER A1\r')
##        self.s.send('SER A2\r')
##        self.s.send('SER A3\r')
##        self.s.send('SER A4\r')
##        self.s.send('SER A5\r')
        time.sleep(1.)
        print self.s.recv(1024)
        

    def recv(self):
        self.s.recv(1024)
        
    def close(self):
        self.s.close()
        
    def move(self,stage,steps):
        if steps < 5000:
            self.s.send('REL A'+str(stage)+' '+str(steps)+' G\r')
            print self.s.recv(1024)
        else:
            print 'Too many steps...'
            
    def pos(self):
        self.s.send('POS\r')
        print self.s.recv(1024)

    def custom(self,command):
        self.s.send(command)
        print self.s.recv(1024)

    def stop(self):
        self.s.send('STO\r')
        print self.s.recv(1024)

    def res(self,fine=True):
        if fine is True:
            self.s.send('RES FINE\r')
        else:
            self.s.send('RES COARSE\r')
        print self.s.recv(1024)

class PSPC:
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def position(x,y):
        print 'X: ' - str(self.x+10.87*x)
        print 'Y: ' + str(self.y+10.87*y)

    def setorigin(x,y):
        self.x=x
        self.y=y

class PIXI:
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def position(x,y):
        print x
        ###print 'X: ' - str(self.x+10.87*x)
        ###print 'Y: ' + str(self.y+10.87*y)
    def setorigin(x,y):
        self.x=x
        self.y=y    

#From three points, fit a circle
#Assume all points are positive
def circ(x,rad,x0,y0):
    return sqrt(rad**2 - (x-x0)**2) + y0
def fitradius(x,y,rad,x0,y0):
    fit = curve_fit(circ,x,y,p0=[rad,x0,y0])
    return fit

def findcentroid(d):
    x,y = meshgrid(arange(shape(d)[1]),arange(shape(d)[0]))
    cx = nansum(x*d)/nansum(d)
    cy = nansum(y*d)/nansum(d)
    return cy,cx

def findcoords(d,bx1,bx2,by1,by2):
    cx,cy = findcentroid(d[bx1:bx2,by1:by2])
    cx = (bx1+cx)*10.84e-3*8
    cy = (by1+cy)*10.84e-3*8
    return cx,cy

#Double Gaussian plus background
def func_twogauss(x, norm1, cent1, width1, norm2, cent2, width2, cons):
    return norm1*exp(-(x-cent1)**2/(2*width1**2)) + cons + \
           norm2*exp(-(x-cent2)**2/(2*width2**2))
