# -*- coding: utf-8 -*-
import numpy as np
import pyfits

nchan = 1024
chan = np.arange(nchan)
mcahist = 0.0*chan
acqaccum = 300.0

# first column is the array of channel numbers 
col1 = pyfits.Column(name='CHANNEL', format='I', array=chan)
# second column is the array of counts per channel
col2 = pyfits.Column(name='COUNTS', format='E', array=mcahist)
# combine the columns to make a table header data unit
tbhdu = pyfits.new_table(pyfits.ColDefs([col1, col2]))
# add keywords to the table header
tbhdu.header.update('DETCHANS', nchan, 'Total number of detector channels available')
tbhdu.header.update('TLMIN1', 0, 'Lowest legal channel number')
tbhdu.header.update('TLMAX1', nchan-1, 'Highest legal channel number')
tbhdu.header.update('EXPOSURE', acqaccum, 'exposure/accumulation time')
tbhdu.header.update('EXTNAME', 'SPECTRUM', 'name of this binary table extension')
tbhdu.header.update('TELESCOP', 'GEMS', 'Telescope (mission) name')
tbhdu.header.update('INSTRUME', 'MOXTEK-XPIN/MCAPIC', 'Instrument name')
tbhdu.header.update('CHANTYPE', 'PHA', 'Channels assigned by detector electronics')
# add DATE-OBS should be in UT-TIME, DATE-END
# 2006-01-13T19:03:57 Start time (UTC) of exposure
# time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime(t))  t = time.time()
# TIMESYS TT
# CREATOR, name of program
# DATE File creation date, same format as DATE-OBS 
# write as a file
tbhdu.writeto('table.fits')


# if multiple tables are wanted, need to define header
# and add each table
# define the main header
#hdu = pyfits.PrimaryHDU(n)
#thdulist = pyfits.HDUList([hdu, tbhdu])
#thdulist.writeto(’table.fits’)

