from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.pyplot import *

import Tkinter as Tk
import sys

import u3
import time
from numpy import *
import pdb
import threading

def convsec(inp):
    hours = inp/60.**2
    minutes = (inp % 60.**2) / 60.
    sec = inp % 60.
    res = '%i:%02i:%02i' % (trunc(hours),trunc(minutes),sec)
    return res
    
class readLabJack(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
    def run(self):
        global starttime
        global saved
        global stopflag

        global tempdata
        global timedata
        global presdata

        #Open up LabJack
        d = u3.U3()
        d.configIO(FIOAnalog=0xff)

        timedata = array(0) #Initialize time array
        tempdata = array((d.getAIN(2))*100-273.15) #Initialize temp array
        presdata = d.getAIN(0)*1000. #Initialize pressure array

        starttime = time.time()
        stime = time.time()
        #Test loop
        while stopflag==0:
            ctime = time.time() #Update time counter
            timedata = append(timedata,ctime-starttime) #Add current elapsed
                                                        #time to data array
            temp = ((d.getAIN(2))*100-273.15) #Update temp
            tempdata = append(tempdata,temp)  #Add current temp to data array
            #pressure = 10.**(d.getAIN(0)-4.) #Update pressure
            pressure = d.getAIN(0)*1000.
            presdata = append(presdata,pressure) #Add current pressure
                                                 #to data array
            time.sleep(1)
            if (time.time()-stime) > 3600.:
                ltime = time.localtime()
                savetxt('ThermalTest_'+str(ltime[0])+'_'+str(ltime[1])+'_'+str(ltime[2])+'_'\
                    +str(ltime[3])+'_'+str(ltime[4])+'.txt',\
                    transpose([timedata,tempdata,presdata]))
                _plot()
                stime = time.time()
                timedata = timedata[-1]
                tempdata = tempdata[-1]
                presdata = presdata[-1]

        #Close LabJack and save data
        d.close()
        savetxt('ThermalTestResults.txt',transpose([timedata,tempdata,presdata]))
        saved = 1

def _plot():
    while size(timedata)!=size(tempdata):
        pass
    tempplt.plot(timedata/(60.**2),tempdata,'r.')
    while size(timedata)!=size(presdata):
        pass
    presplt.plot(timedata/(60.**2),presdata,'r.')
    canvas.show()

def _quit():
    global stopflag
    stopflag = 1

def _start():
    global stopflag
    global starttime
    global saved
    
    global tempdata
    global presdata
    global timedata

    #Disable button
    start.config(state='disabled')

    #Start LabJack acquisition
    read1 = readLabJack()
    read1.start()
    saved = 0

    #Enter loop to update labels and detect outside events
    time.sleep(2) #Wait for acquisition to populate vectors
    while stopflag == 0:
        time.sleep(1)
        tempstr.set('Temperature: %.2f' % tempdata[-1]) #Update temp label
        presstr.set('Pressure: %.2f' % presdata[-1]) #Update pressure label
        timestr.set('Elapsed Time: ' + convsec(timedata[-1])) #Update time label
        root.update() #Detect outside events

    #Wait for data to be saved before exiting
    while saved==0:
        pass
        
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate

root = Tk.Tk()
root.wm_title("Window Thermal Test")

#Sets up subplots
#matplotlib.rcParams['axes.labelsize']=2
f = Figure(figsize=(7,6))
f.hold(True)
tempplt = f.add_subplot(211)
presplt = f.add_subplot(212)
setp(tempplt.get_xticklabels(),visible=False)
tempplt.set_ylabel('Temperature (C)')
presplt.set_ylabel('Pressure (Torr)')
presplt.set_xlabel('Elapsed Time (hrs)')
tempplt.set_title('LabJack Readout')

global stopflag
stopflag=0

# a tk.DrawingArea
canvas = FigureCanvasTkAgg(f, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=0)

#toolbar = NavigationToolbar2TkAgg( canvas, root )
#toolbar.update()
canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=0)

#Interactive frame
ifrm = Tk.Frame(master=root)
ifrm.pack(side=Tk.BOTTOM,fill=Tk.BOTH,expand=1)
#Frame for Labels
lfrm = Tk.Frame(master=ifrm)
lfrm.pack(side=Tk.LEFT,fill=Tk.X,expand=0)
#Frame for Buttons
bfrm = Tk.Frame(master=ifrm)
bfrm.pack(side=Tk.LEFT,fill=Tk.Y,expand=0)

#Buttons
stop = Tk.Button(master=bfrm,text='Stop',command=_quit)
stop.pack(side=Tk.BOTTOM)
start = Tk.Button(master=bfrm,text='Start',command=_start)
start.pack(side=Tk.BOTTOM)
visualize = Tk.Button(master=bfrm,text='Plot',command=_plot)
visualize.pack(side=Tk.RIGHT)

#Info Labels
presstr = Tk.StringVar()
preslabel = Tk.Label(master=lfrm,textvariable=presstr,width=70)
preslabel.pack(side=Tk.TOP)
tempstr = Tk.StringVar()
templabel = Tk.Label(master=lfrm,textvariable=tempstr,width=70)
templabel.pack(side=Tk.TOP)
timestr = Tk.StringVar()
timelabel = Tk.Label(master=lfrm,textvariable=timestr,width=70)
timelabel.pack(side=Tk.TOP)

Tk.mainloop()
