import numpy as np
import matplotlib.pyplot as plt
import utilities.imaging.analysis as anal
import os,glob,pyfits

os.chdir('/home/rallured/Dropbox/Arcus/OGRETest/')

def getdata():
    fn = glob.glob('*.fits')
