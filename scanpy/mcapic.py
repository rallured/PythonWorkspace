""" code to run MCAPIC
In ipython load using "from mcapic import *"
then init() to setup, acquire(time) to acquire, close() to end
"""

import serial, time, struct, pickle, pyfits
import numpy as np
import matplotlib.pyplot as plt

def plotmca(roi=[]):
  plt.clf()
  plt.plot(chan, mcahist, 'k-') # plot as black line
  ymax = 1.1*max(max(mcahist), 10)
  plt.axis([0,nchan,0, ymax])
  plt.xlabel('Channel')
  plt.ylabel('Number events')
  plt.title('MCA Data')
  if (len(roi) == 2):
    plt.plot([roi[0], roi[0]], [0, ymax])
    plt.plot([roi[1], roi[1]], [0, ymax])
  plt.draw()

def init():
  global nchan, chan, mcahist, ser, plt
  nchan = 512# number of MCA channels
  chan = np.arange(nchan) # array of channel numbers
  mcahist = np.zeros(nchan) # array to store MCA data
  # open serial communications with MCAPIC
  ser = serial.Serial('COM1', 115200, timeout=1)
  ser.timeout = 1 # set one second time out for serial input from mcapic
  # set up plot window for MCA data
  plt.ion()
  plotmca()
  print "DAQ set up"

def acquire(acqtime, roi=[]):
  global nchan, mcahist, ser, plt, realtime, acqaccum, starttime, stoptime
  print "Starting acquisition",
  for i in range(nchan):
      mcahist[i] = 0
  acqinc = 10.0 # acquistion time increment (seconds)
  acqaccum = 0.0 # accumulated acquisition time (seconds)
  starttime = time.time()
  plotmca()
  if sum(mcahist) > 0:
    print "mcahist not zeroed"
    mcahist[:] = 0.0
  # remove any previous data
  #ser.flushInput() # discard anything sent previously by mcapic
  #ser.write('E') # stop any previous acquisition
  #a = ser.read(2*nchan) # read the data, and then discard
  while acqaccum < acqtime:
    ser.write('S') # start the MCAPIC
    time.sleep(acqinc) # wait
    ser.write('E') # stop the acquisition
    a = ser.read(2*nchan) # read the data
    if len(a) <> 2*nchan:
        print "Did not read correct number of bytes from mcapic ", len(a)
    b = struct.unpack("1024B",a) # unpack the data from string to bytes
    c = np.zeros(nchan)
    for i in range(nchan):
        c[i] = b[2*i]+256*b[2*i+1]

###rdaly::added code to test new timer functionality
    timer_count = ser.read(4) # read the timer data
    lsw_lower,lsw_upper,msw_lower,msw_upper = struct.unpack("4B",timer_count)

    seconds_lsw = (lsw_upper*256 + lsw_lower)
    seconds_msw = (msw_upper*256 + msw_lower)
    seconds_acq = (seconds_msw*65536 + seconds_lsw) / 156250.0
###rdaly::end code added to test new timer functionality

    print "sum(c)=", sum(c), 
    #print c

###rdaly::added print statement for timer diagnostic
    print " ", seconds_acq, " seconds"
###rdaly::end code added for diagnostic

    mcahist = mcahist+c
    #for i in range(nchan):
    #    mcahist[i] += b[2*i]+256*b[2*i+1]
    acqaccum += seconds_acq
    # print statistics
    print "\rTime=%d  events=%d" % (acqaccum, sum(mcahist)),
    # print statistics for any regions of interest
    if (len(roi) == 2):
      q = (chan >= roi[0]) & (chan <= roi[1])
      s = sum(mcahist[q])
      sx = sum(mcahist[q]*chan[q])
      sx2 = sum(mcahist[q]*chan[q]*chan[q])
      print " ROI events=%d mean=%f" % (s, sx/s),
      plotmca(roi)
    else :
      plotmca() # plot the data
  stoptime = time.time()-starttime
  print

def save_pickle(filename):
  global nchan, mcahist, realtime, acqaccum
  # save the data  
  mcafile = open(filename, 'wb')
  pickle.dump(acqaccum, mcafile)
  pickle.dump(realtime, mcafile)
  pickle.dump(nchan, mcafile)
  pickle.dump(mcahist, mcafile)
  mcafile.close()

def save(filename, keywords=[]):
# save the data in a FITS file
# write any keywords in the lis keywords to the header
  global nchan, chan, mcahist, realtime, acqaccum, starttime, stoptime
  # first column is the array of channel numbers 
  col1 = pyfits.Column(name='CHANNEL', format='E', array=chan)
  # second column is the array of counts per channel
  col2 = pyfits.Column(name='COUNTS', format='E', array=mcahist)
  # combine the columns to make a table header data unit
  tbhdu = pyfits.new_table(pyfits.ColDefs([col1, col2]))
  # add any user specified keywords to the header
  for i in range(len(keywords)):
    t = keywords[i]
    tbhdu.header.update(t[0], t[1], t[2])
  # add keywords to the table header
  tbhdu.header.update('DETCHANS', nchan, 'Total number of detector channels available')
  tbhdu.header.update('TLMIN1', 0, 'Lowest legal channel number')
  tbhdu.header.update('TLMAX1', nchan-1, 'Highest legal channel number')
  tbhdu.header.update('EXPOSURE', acqaccum, 'exposure/accumulation time')
  tbhdu.header.update('EXTNAME', 'SPECTRUM', 'name of this binary table extension')
  tbhdu.header.update('TELESCOP', 'GEMS', 'Telescope (mission) name')
  tbhdu.header.update('INSTRUME', 'MOXTEK-XPIN/MCAPIC', 'Instrument name')
  tbhdu.header.update('CHANTYPE', 'PHA', 'Channels assigned by detector electronics')
  t = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime(starttime))
  tbhdu.header.update('DATE-OBS', t, 'Start time (UTC) of exposure')
  t = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime(stoptime))
  tbhdu.header.update('DATE-END', t, 'End time (UTC) of exposure')
  t = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime(time.time()))
  tbhdu.header.update('DATE', t, 'File creation date (UTC)')
  tbhdu.header.update('TIMESYS', 'TT')
  tbhdu.header.update('CREATOR', 'mcapic.py (ver. 2010-03-29)')
  tbhdu.writeto(filename)

def close():
  global ser, plt
  ser.close() # close the serial port
  plt.close() # close the plot window
  print "DAQ closed"

