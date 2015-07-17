'''calcola riflettivita' per oggetto wolter con
geometria ottenuta dal file imp_offAxis col nome passato per
argomento. (per essere chiamato da VB)'''

import sys
from userKov.pyGeoShell.CsimpleCreaDiam import CSimpleCreaDiam
from pylab import *
#import dislin
from userKov.pyML.Ccoating import Ccoating        
    

    
def wolterTel(geoFile="imp_OffAxis.txt",referFile=0):
    '''calcola area efficace totale per telescopio wolter da file di
    impostazioni della geometria'''
    #tel=CSimpleCreaDiam (sys.argv[1])
    a=CSimpleCreaDiam(geoFile)
    ener=map(float,range(1,51))
    area=a.aeff(ener)
    '''if referFile:
        refer=load(referFile)
        referX,referY=refer[0,:],refer[1,:]
        plot(ener,area,"-",referX,referY,"-")
    else:
        plot(ener,area,"-")
    show()'''
    
    open("diamWolter.dat","w").write("\n".join([str(d) for d in a.d]))
    open("results.txt","w").write("\n".join(["%s\t%s" %(ee,aa) for ee,aa in zip(ener, area)]))
    return a 

def testShell(angle,ener):
    '''plotta e scrive su file la 
    riflettività di un coating a un dato angolo e vettore energie'''
    shell=Ccoating() #inizializza coating con ricetta di default
    shell.ener=ener
    reflex=shell.reflex(angle)
    #plot(ener,reflex)
    #show()
##    dislin.plot(ener,reflex)
##    dislin.disfin()
    open("reflexAt%03i.txt"%(int(angle*1000)),"w").write("\n".join(["%s\t%s" %(ee,aa) for ee,aa in zip(ener, reflex)]))


if __name__=="__main__":

    #ener=[i/10. for i in range(1,200)]
    #angle=0.006
    #testShell(angle,ener)
    referFile="sommaTOT_sxbl.txt"
    w=wolterTel("imp_OffAxis_blIR.txt",referFile)

