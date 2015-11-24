import numpy as np
import matplotlib.pyplot as plt
import pyfits,os


#Load in IRIS and TFT IFs
os.chdir('/home/rallured/Dropbox/WFS/SystemAlignment/TFT1/')
iris = pyfits.getdata('IRIS2/151014_IRIS_IFs.fits')
tft = pyfits.getdata('TFTActuation/151014_TFT_IFs.fits')

def plotIF(i):
    """Plot IRIS and TFT IFs side by side"""
    fig = plt.figure(figsize=(15,6))
    fig.add_subplot(121)
    plt.imshow(iris[i])
    plt.colorbar()
    plt.title('IRIS Cell '+str(i+1))
    fig.add_subplot(122)
    plt.imshow(tft[i])
    plt.colorbar()
    plt.title('TFT Cell '+str(i+1))
    n = i+1
    plt.savefig('Cell%02i.png'%n)
