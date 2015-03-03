#This submodule contains image analysis routines
from matplotlib.pyplot import *
from numpy import *

def findMoments(d):
    x,y = meshgrid(arange(shape(d)[1]),arange(shape(d)[0]))
    cx = nansum(x*d)/nansum(d)
    cy = nansum(y*d)/nansum(d)
    rmsx = nansum((x-cx)**2*d)/nansum(d)
    rmsy = nansum((y-cy)**2*d)/nansum(d)
    
    return cx,cy,sqrt(rmsx),sqrt(rmsy)

class pointGetter:
    """Creates an object tied to an imshow where the user can
    accumulate a list of x,y coords by right clicking on them.
    When done, the user can press the space key.
    """
    def __init__(self,img):
        self.fig = figure()
        self.x = zeros(0)
        self.y = zeros(0)
        self.con = self.fig.canvas.mpl_connect('key_press_event',\
                                                self.keyEvent)
        ax = self.fig.add_subplot(111)
        ax.imshow(img)
    def keyEvent(self,event):
        if event.key==' ':
            self.x = append(self.x,event.xdata)
            self.y = append(self.y,event.ydata)
    def close(self):
        self.fig.canvas.mpl_disconnect(self.con)
        close(self.fig)

class peakInformation:
    def __init__(self,img):
        self.img = img
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.gx = []
        self.gy = []
        self.rmsx = []
        self.rmsy = []
        self.fig = figure()
        ax = self.fig.add_subplot(111)
        ax.imshow(img)
        self.con = self.fig.canvas.mpl_connect('button_press_event',self.clickEvent)
    def clickEvent(self,event):
        #If x0 and y0 are undefined, set them and return
        if self.x0 is None:
            #Define first point
            self.x0 = event.xdata
            self.y0 = event.ydata
            return
        #If x1 and y1 are undefined, set them and return centroid
        #Define second point
        self.x1 = event.xdata
        self.y1 = event.ydata
        #Order points properly
        x0 = min([self.x0,self.x1])
        x1 = max([self.x0,self.x1])
        y0 = min([self.y0,self.y1])
        y1 = max([self.y0,self.y1])            
        #Compute centroid between coordinates and update centroid list
        cx,cy,rmsx,rmsy = findMoments(self.img[y0:y1,x0:x1])
        print 'X: ' + str(cx+x0)
        print 'Y: ' + str(cy+y0)
        print 'RMS X: ' + str(rmsx)
        print 'RMS Y: ' + str(rmsy)
        try:
            self.gx.append(cx+x0)
            self.gy.append(cy+y0)
            self.rmsx.append(rmsx)
            self.rmsy.append(rmsy)
        except:
            self.gx = [cx+x0]
            self.gy = [cy+y0]
            self.rmsx = [rmsx]
            self.rmsy = [rmsy]
        #Reset event
        self.x0 = None
        self.x1 = None
        self.y0 = None
        self.y1 = None
    def close(self):
        self.fig.canvas.mpl_disconnect(self.con)
        print self.gx
        print self.gy
        print self.rmsx
        print self.rmsy
        close(self.fig)

def getPoints(img):
    """This function will bring up the image, wait for the user to
    right click on any number of points, then wait for the user to
    end with the 'q' key, and then return a list of x,y coords
    """
    #Create instance of pointGetter class
    p = pointGetter(img)

    #When user presses enter, close the pointGetter class and
    #return the list of coordinates
    try:
        input("Press enter when finished selecting points...")
    except:
        pass
    x = p.x
    y = p.y
    p.close()

    return x,y
