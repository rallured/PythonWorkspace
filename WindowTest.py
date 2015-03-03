import u3, u6
import time
from datetime import datetime
import pdb
from numpy import *
import os
from Tkinter import *
import tkMessageBox
import sys
import thread

def Scan():
    global temp
    global pressure
    global etime0
    global etime1
    #Open LabJack and set FIO4-7 to analog
    d = u3.U3()
    d.configIO(FIOAnalog=15)

    #Configure stream
    d.streamConfig(NumChannels=4,PChannels=[0,2,200,224],NChannels=[31,31,31,31],\
                   Resolution=0,SampleFrequency=10)
    numrecords = 0

    #Start stream
    starttime = time.time()
    d.streamStart()

    #Read in data
    for r in d.streamData():
        if numrecords > 20:
            break
        try:
            temp.extend(r['AIN2'])
            pressure.extend(r['AIN0'])
            etime0.extend(r['AIN200'])
            etime1.extend(r['AIN224'])
            numrecords = size(temp)
        except:
            temp = r['AIN2']
            pressure = r['AIN0']
            etime0 = r['AIN200']
            etime1 = r['AIN224']

    print time.time()-starttime
    d.streamStop()

    #Close LabJack
    d.close()

    #Save results
    os.chdir('/Users/ryanallured/PythonWorkspace/Testing/')
    savetxt('ScanResults.txt',[temp,pressure])

thread.start_new_thread(Scan,())
print 'test'
