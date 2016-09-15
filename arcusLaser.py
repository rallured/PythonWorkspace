import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import tkFileDialog
import threading

#Global variables
drefcen = (0,0)
dcurcen = (0,0)
rrefcen = (0,0)
rcurcen = (0,0)

#Calibration factors
calP = 0.
calR = 0.
calY = 0.

#Reference images
rrefimg = None
drefimg = None
rcurimg = None
dcurimg = None

def yaw():
    return (dcurcen[0]-drefcen[0])*calY

def pitch():
    return (rcurcen[0]-rrefcen[0])*calP

def roll():
    return (rcurcen[1]-rrefcen[1])*calR

def getDiffraction(f=None):
    """
    Dummy function that gets image from diffraction camera.
    Second argument returns coordinates of centroid
    """
    img = np.random.random(size=(128,128))
    return img,(100,60)

def getReflection(f=None):
    """
    Dummy function that gets image from reflection camera.
    Second argument returns coordinates of centroid
    """
    img = np.random.random(size=(128,128))
    return img,(100,80)

def _dref():
    """
    Acquire the diffraction reference image and centroid
    """
    #Update plot
    img,cen = getDiffraction()
    dref.imshow(img)
    dref.plot(cen[0],cen[1],'*b',markersize=14)
    dref.set_ylim([0,128])
    dref.set_xlim([0,128])
    canvas.show()

    #Set global variables
    global drefcen,drefimg
    drefcen = cen
    drefimg = img

    return

def _drefF():
    """
    Acquire the diffraction reference image and centroid FROM FILE
    """
    #Get filename
    f = tkFileDialog.askopenfilename()
    
    #Update plot
    img,cen = getDiffraction(f=f)
    dref.imshow(img)
    dref.plot(cen[0],cen[1],'*b',markersize=14)
    dref.set_ylim([0,128])
    dref.set_xlim([0,128])
    canvas.show()

    #Set global variables
    global drefcen,drefimg
    drefcen = cen
    drefimg = img
    
    return

def _drefsave():
    """Save the diffraction reference to file"""
    #Get filename
    f = tkFileDialog.asksaveasfilename()
    #Save image
    try:
        pyfits.writeto(f,drefimg)
    except:
        _drefsave()
    return

def _dcur():
    """
    Acquire the diffraction current image and centroid
    Calculate yaw and update label
    """
    #Update plot
    img,cen = getDiffraction()
    dcur.imshow(img)
    dcur.plot(cen[0],cen[1],'*b',markersize=14)
    dcur.set_ylim([0,128])
    dcur.set_xlim([0,128])
    canvas.show()

    #Set global variables
    global dcurcen,dcurimg
    dcurcen = cen
    dcurimg = img

    #Update yaw label
    yawstr.set('Yaw (Diffraction): '+str(round(yaw())))

    return

def _dcurF():
    """
    Acquire the diffraction current image and centroid FROM FILE
    """
    #Get filename
    f = tkFileDialog.askopenfilename()
    
    #Update plot
    img,cen = getDiffraction(f=f)
    dcur.imshow(img)
    dcur.plot(cen[0],cen[1],'*b',markersize=14)
    dcur.set_ylim([0,128])
    dcur.set_xlim([0,128])
    canvas.show()

    #Set global variables
    global dcurcen,dcurimg
    dcurcen = cen
    dcurimg = img

    #Update yaw label
    yawstr.set('Yaw (Diffraction): '+str(round(yaw())))
    
    return

def _dcursave():
    """Save the diffraction current to file"""
    #Get filename
    f = tkFileDialog.asksaveasfilename()
    #Save image
    try:
        pyfits.writeto(f,dcurimg)
    except:
        _dcursave()
    return

def _rref():
    """
    Acquire the reflection reference image and centroid
    """
    #Update plot
    img,cen = getReflection()
    rref.imshow(img)
    rref.plot(cen[0],cen[1],'*b',markersize=14)
    rref.set_ylim([0,128])
    rref.set_xlim([0,128])
    canvas.show()

    #Set global variables
    global rrefcen,rrefimg
    rrefcen = cen
    rrefimg = img

    return

def _rrefF():
    """
    Acquire the reflection reference image and centroid FROM FILE
    """
    #Get filename
    f = tkFileDialog.askopenfilename()
    #Update plot
    img,cen = getReflection(f=f)
    rref.imshow(img)
    rref.plot(cen[0],cen[1],'*b',markersize=14)
    rref.set_ylim([0,128])
    rref.set_xlim([0,128])
    canvas.show()

    #Set global variables
    global rrefcen,rrefimg
    rrefcen = cen
    rrefimg = img

    return

def _rrefsave():
    """Save the reflection reference to file"""
    #Get filename
    f = tkFileDialog.asksaveasfilename()
    #Save image
    try:
        pyfits.writeto(f,rrefimg)
    except:
        _rrefsave()
    return

def _rcur():
    """
    Acquire the diffraction current image and centroid
    Calculate yaw and update label
    """
    #Update plot
    img,cen = getReflection()
    rcur.imshow(img)
    rcur.plot(cen[0],cen[1],'*b',markersize=14)
    rcur.set_ylim([0,128])
    rcur.set_xlim([0,128])
    canvas.show()

    #Set global variables
    global rcurcen,rcurimg
    rcurcen = cen
    rcurimg = img

    #Update labels
    pitchstr.set('Pitch (Reflection): '+str(round(pitch())))
    rollstr.set('Roll (Reflection): '+str(round(roll())))

    return

def _rcurF():
    """
    Acquire the reflection current image and centroid FROM FILE
    """
    #Get filename
    f = tkFileDialog.askopenfilename()
    #Update plot
    img,cen = getReflection(f=f)
    rcur.imshow(img)
    rcur.plot(cen[0],cen[1],'*b',markersize=14)
    rcur.set_ylim([0,128])
    rcur.set_xlim([0,128])
    canvas.show()

    #Set global variables
    global rcurcen,rcurimg
    rcurcen = cen
    rcurimg = img

    #Update labels
    pitchstr.set('Pitch (Reflection): '+str(round(pitch())))
    rollstr.set('Roll (Reflection): '+str(round(roll())))

    return

def _rcursave():
    """Save the reflection current to file"""
    #Get filename
    f = tkFileDialog.asksaveasfilename()
    #Save image
    try:
        pyfits.writeto(f,rcurimg)
    except:
        _rcursave()
    return
    
#Set up Tk
root = Tk.Tk()
root.wm_title('Arcus Yaw Laser')

#Image figure
imgs = plt.figure(figsize=(12,10))
imgs.hold(True)
dref = imgs.add_subplot(221)
dref.set_title('Diffraction Reference')
dcur = imgs.add_subplot(222)
dcur.set_title('Diffraction Current')
rref = imgs.add_subplot(223)
rref.set_title('Reflection Reference')
rcur = imgs.add_subplot(224)
rcur.set_title('Reflection Current')

# a tk.DrawingArea
canvas = FigureCanvasTkAgg(imgs, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)
canvas.show()

#Frame for information
infofrm = Tk.Frame(master=root)
infofrm.pack(side=Tk.BOTTOM,fill=Tk.BOTH,expand=1)

#Frame for reference labels
reffrm = Tk.Frame(master=infofrm)
reffrm.pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)
reffrm2 = Tk.Frame(master=infofrm)
reffrm2.pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)

#Frame for current labels
curfrm = Tk.Frame(master=infofrm)
curfrm.pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)
curfrm2 = Tk.Frame(master=infofrm)
curfrm2.pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)

#Frame for angles
angfrm = Tk.Frame(master=infofrm)
angfrm.pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)

#Info labels
drefstr = Tk.StringVar()
drefstr.set('Diffraction Ref Centroid: ')
dreflbl = Tk.Label(master=reffrm,textvariable=drefstr,width=40)
dreflbl.config(font=('Courier',20),justify='left')
dreflbl.pack(side=Tk.TOP)

dcurstr = Tk.StringVar()
dcurstr.set('Diffraction Cur Centroid: ')
dcurlbl = Tk.Label(master=curfrm,textvariable=dcurstr,width=40)
dcurlbl.config(font=('Courier',20),justify='left')
dcurlbl.pack(side=Tk.TOP)

rrefstr = Tk.StringVar()
rrefstr.set('Reflection Ref Centroid: ')
rreflbl = Tk.Label(master=reffrm,textvariable=rrefstr,width=40)
rreflbl.config(font=('Courier',20),justify='left')
rreflbl.pack(side=Tk.TOP)

rcurstr = Tk.StringVar()
rcurstr.set('Reflection Cur Centroid: ')
rcurlbl = Tk.Label(master=curfrm,textvariable=rcurstr,width=40)
rcurlbl.config(font=('Courier',20),justify='left')
rcurlbl.pack(side=Tk.TOP)

pitchstr = Tk.StringVar()
pitchstr.set('Pitch (Reflection): ')
pitchlbl = Tk.Label(master=angfrm,textvariable=pitchstr,width=40)
pitchlbl.config(font=('Courier',20),justify='left')
pitchlbl.pack(side=Tk.TOP)

rollstr = Tk.StringVar()
rollstr.set('Roll (Reflection): ')
rolllbl = Tk.Label(master=angfrm,textvariable=rollstr,width=40)
rolllbl.config(font=('Courier',20),justify='left')
rolllbl.pack(side=Tk.TOP)

yawstr = Tk.StringVar()
yawstr.set('Yaw (Diffraction): ')
yawlbl = Tk.Label(master=angfrm,textvariable=yawstr,width=30)
yawlbl.config(font=('Courier',20))
yawlbl.pack(side=Tk.TOP)

calPstr = Tk.StringVar()
calPstr.set('Pitch Calibration: %.2f' % calP)
calPlbl = Tk.Label(master=angfrm,textvariable=calPstr,width=40)
calPlbl.config(font=('Courier',20),justify='left')
calPlbl.pack(side=Tk.TOP)

calRstr = Tk.StringVar()
calRstr.set('Roll Calibration: %.2f' % calR)
calRlbl = Tk.Label(master=angfrm,textvariable=calRstr,width=40)
calRlbl.config(font=('Courier',20),justify='left')
calRlbl.pack(side=Tk.TOP)

calYstr = Tk.StringVar()
calYstr.set('Yaw Calibration: %.2f' % calY)
calYlbl = Tk.Label(master=angfrm,textvariable=calYstr,width=30)
calYlbl.config(font=('Courier',20))
calYlbl.pack(side=Tk.TOP)

#Buttons
drefb = Tk.Button(master=reffrm,text='D Ref',command=_dref)
drefb.pack(side=Tk.LEFT)
drefbf = Tk.Button(master=reffrm2,text='D Load File',command=_drefF)
drefbf.pack(side=Tk.LEFT)
drefbs = Tk.Button(master=reffrm2,text='D Save File',command=_drefsave)
drefbs.pack(side=Tk.LEFT)
rrefb = Tk.Button(master=reffrm,text='R Ref',command=_rref)
rrefb.pack(side=Tk.LEFT)
rrefbf = Tk.Button(master=reffrm2,text='R Load File',command=_rrefF)
rrefbf.pack(side=Tk.LEFT)
rrefbs = Tk.Button(master=reffrm2,text='R Save File',command=_rrefsave)
rrefbs.pack(side=Tk.LEFT)
dcurb = Tk.Button(master=curfrm,text='D Cur',command=_dcur)
dcurb.pack(side=Tk.LEFT)
dcurbf = Tk.Button(master=curfrm2,text='D Load File',command=_dcurF)
dcurbf.pack(side=Tk.LEFT)
dcurbs = Tk.Button(master=curfrm2,text='D Save File',command=_dcursave)
dcurbs.pack(side=Tk.LEFT)
rcurb = Tk.Button(master=curfrm,text='R Cur',command=_rcur)
rcurb.pack(side=Tk.LEFT)
rcurbf = Tk.Button(master=curfrm2,text='R Load File',command=_rcurF)
rcurbf.pack(side=Tk.LEFT)
rcurbs = Tk.Button(master=curfrm2,text='R Save File',command=_rcursave)
rcurbs.pack(side=Tk.LEFT)

root.mainloop()
