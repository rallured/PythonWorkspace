import serial, time, struct, pickle

nchan = 1024 # number of channels
acqtime = 60 # acquistion time (seconds)

ser = serial.Serial('COM1', 115200, timeout=1)
print "Start acquistion and wait ", acqtime, " seconds"
ser.write('S') # start the MCAPIC
time.sleep(acqtime) # wait
ser.write('E') # stop the acquisition
print "Acquistion done"
a = ser.read(2*nchan) # read the data
print "Data read"
ser.close() # close the serial port
print "Serial port closed"
b = struct.unpack("2048B",a)
mcaout = nchan*[None] # array to store histogram
for i in range(nchan):
  mcaout[i] = b[2*i]+256*b[2*i+1]

print mcaout
mcafile = open('mca.dat', 'wb')
pickle.dump(mcaout, mcafile)
mcafile.close()
