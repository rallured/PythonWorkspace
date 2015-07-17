from generalRoutines import interpola
from generalRoutines import loadCol
import os
from pylab import *

def applyFilter(xFilter,yFilter,xData,yData):
    s=[]
    if not iterable(yData):
        xData=[xData]
        yData=[yData]         
    for xD,yD in zip(xData,yData):
        yInt=interpola (xFilter,yFilter,xD)
        yD=array(yInt)*yD
        for e1,e2 in zip(xD,yD):
            s.append("%s\t%s"%(e1,e2))
        s.append("\n")
    return "\n".join(s)

def filterOaAree(dataFile,filterFile,outFile=None):
    print "filtro %s"%dataFile
    xFilter=load(filterFile)[:,0]	
    yFilter=load(filterFile)[:,1]
    xData=loadCol(dataFile,0)
    yData=loadCol(dataFile,1)
    xFilter=xFilter/1000
    toWrite=applyFilter(xFilter,yFilter,xData,yData)
    #plot (x,y,'-')
    #show()
    if not outFile:outFile="_QE".join(os.path.splitext(dataFile))
    f=open(outFile,"w")
    f.write(toWrite)
    f.close()        
        

if __name__=="__main__":
    
    folder=r"C:\work\copiaOA\edge_eq3\thermalCoating\spider_filter_qe_ufficiali"
    filt=r"C:\work\copiaOA\edge_eq3\thermalCoating\spider_filter_qe_ufficiali\mos2_qe.qdp"
    listaFile=r"xQE.txt"
    for dataFile in open(os.path.join(folder,listaFile),"r").read().split():
        filterOaAree(os.path.join(folder,dataFile),filt)


    