'''legge psfFile generato da OffAxis, fa conti'''

import matplotlib
from pylab import *
import numpy
import extractInfo
import os

def hystogram(ar):
    selectedDic={}
    for ang1 in ar:
        if selectedDic.has_key(ang1): selectedDic[ang1]=selectedDic[ang1]+1
        else: selectedDic[ang1]=1
    x=selectedDic.keys()										
    x.sort()
    y=[selectedDic[xx] for xx in x]
    return x,y

def plotDist(fpFile):     
    t=load (fpFile)
    selected=compress(t[:,2]==1.,t,0) #prende gli elementi con qa=1 (doppia riflessione)
    selected=compress(selected[:,8]!=0,selected,0) #prende gli elementi che non hanno angolo 0 sulla prima shell
    x2,y2=hystogram(selected[:,8])
    #x2,y2=hystogram(selected[:,8])
    #xx2,yy2=hist(selected[:,8])
    #print hist(selected[:,8])
    #plot(x2,y2,"-")   #,xx2,yy2,"-")
    #plot(hist(selected[:,8]))
    #show ()
    return x2,y2

if __name__=="__main__":

    #"test_focalPlane\m2_ff000_195_max65"
    fpFile=r"C:\work\traie4\traie4\provaTraie4\1shell_ref\psf_data_02.txt"  
    plotDist(fpFile)
    
    