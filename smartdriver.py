from numpy import *
from matplotlib.pyplot import *
from ctypes import *
import os

#Define sSDDeviceStruct for interface with Smart Driver DLL
class SDstruct(Structure):
    _fields_ = [
        ("iDriverType",c_int32),
        ("cSerialNumber",c_char*33),
        ("lNIDeviceNumber",c_int32),
        ("lClockSetting",c_int32),
        ("lSDAddress",c_int32),
        ("lNumChannels",c_int32),
        ("fCalibrationArray",(c_long*512)*3),
        ("lCalibrationDate",c_int32)]

#Function to load ctypes float array with values
def fltarrtoc(arr):
    ct = c_float * size(arr) #Define voltage vector array type
    arrc = ct()
    for i in range(size(arr)):
        arrc[i] = arr[i]
    return arrc

class smartDriver:
    def __init__(self):
        #Load DLL
        curdir = os.getcwd()
        os.chdir('C:/Users/WFS/Downloads/Utilities/')
        iris = CDLL('IrisAO')
        os.chdir(curdir)
        #Define pointers for functions
        filename = c_char_p('C:/Users/WFS/Downloads/Utilities/19120001Devices.dcf')
        self.pstruct = pointer(SDstruct())
        self.flag = c_int32(0)
        #Initialize
        init = iris['InitSD']
        init.argtypes = [c_char_p,POINTER(SDstruct),c_int32]
        com = init(filename,self.pstruct,self.flag)
        if com==0:
            print 'Smart Driver initialized...'
        else:
            print com
        #Define set command
        self.setv = iris['SetSDVoltages']
        self.setv.argtypes = [POINTER(SDstruct),POINTER(c_float*128),c_int32]
        #Define shutdown command
        self.shutdown = iris['ShutdownSD']
        self.shutdown.argtypes = [POINTER(SDstruct),c_int32]
        #Initialize voltage array
        self.volt = zeros(128)
    def setvolt(self,voltarr):
        self.volt = voltarr
        com = self.setv(self.pstruct,byref(fltarrtoc(self.volt)),self.flag)
        if com==0:
            print 'Voltages set'
    def setchan(self,chan,voltage):
        self.volt[chan-1] = voltage
        com = self.setv(self.pstruct,byref(fltarrtoc(self.volt)),self.flag)
        if com==0:
            print 'Channel set'
    def shut(self):
        com = self.shutdown(self.pstruct,self.flag)
        if com==0:
            print 'Smart Driver shutdown...'

#Script to check Driver operation
#Sets voltages to 5 volts for 3 seconds
##sd = smartDriver()
##sd.setvolt(repeat(5,128))
##pause(3)
##sd.shut()

#Load DLL
##print os.path.isfile('C:/Users/WFS/Downloads/\
##Iris AO Installation Disk v0.4.0.0 Legacy/\
##Iris AO Installation Disk v0.4.0.0 Legacy/Smart Driver II/Utilities/IrisAO.dll')
##iris = CDLL('C:/Users/WFS/Downloads/\
##Iris AO Installation Disk v0.4.0.0 Legacy/\
##Iris AO Installation Disk v0.4.0.0 Legacy/Smart Driver II/Utilities/IrisAO.dll')
##
###Initialize pointers for DLL functions
##filename = c_char_p('C:/Users/WFS/Downloads/\
##Iris AO Installation Disk v0.4.0.0 Legacy/\
##Iris AO Installation Disk v0.4.0.0 Legacy/Smart Driver II/Utilities/\
##19120001Legacy.dcf')
##pstruct = pointer(SDstruct())
##flag = c_int32(0)
##
###Initialize Smart Driver
##init = iris['InitSD']
##init.argtypes = [c_char_p,POINTER(SDstruct),c_int32]
##com = init(filename,pstruct,flag)
##print com
##
###Create voltage vector - set all channels to 5
##volt = fltarrtoc(repeat(5,128))
##
###Set voltages
##setv = iris['SetSDVoltages']
##setv.argtypes = [POINTER(SDstruct),POINTER(farr),c_int32]
####com = setv(pstruct,volt2.ctypes.data_as(POINTER(c_float*128)),flag)
##com = setv(pstruct,byref(volt),flag)
##print com
##
##pause(3)
##
###Close Smart Driver
##shutdown = iris['ShutdownSD']
##shutdown.argtypes = [POINTER(SDstruct),c_int32]
##com = shutdown(pstruct,flag)
##print com
