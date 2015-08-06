from math import *
from userKov.pyML.Ccoating import Ccoating
from numpy import *
import os
from userKov.pyGeneralRoutines.generalRoutines import plottaFile
from direttoInverso import inverso,nextInner
from userKov.pyGeneralRoutines.generalRoutines import rPars

'''v2.0 26-4-07
modificato calcolo geometria per partire da diametro massimo esterno alla pupilla d'ingresso
con spessore fisso e altezza variabile
vecchie configurazioni di interesse rivalutate con nuovo sistema per verifica.'''

def volume(diam,ang,spes,SH):
    Dmed=diam
    Dmax=Dmed+2*tan(ang)*SH
    Dmin=Dmed-2*tan(3*ang)*SH
    v=((Dmax+Dmin)/2+Dmed)*pi*spes*SH
    return v

def wfxtDiam(Dmin,Dmax,FL,SHpoints,fixedThick,fovArcMin,density,coating=None,ener=None,dHMin=None):
    '''calcola sequenza di diametri. SHpoints sono due coppie di coordinate in grafico (D,SH) nel
    formato ((x1,x2),(y1,y2)). '''
    s=[] #sequenza di stringhe che vengono unite alla fine e restituite come tabella

    fov=radians(float(fovArcMin)/60)
    v=0
    n=1
    totAeff=0
    t=float(fixedThick)    
    SHm,SHq=rPars(SHpoints)
    SH= SHm*Dmax+SHq  #float(SHmin)
    D=inverso(FL,Dmax-2*t,SH)  #float(Dmin)
    if dHMin:
        if D<dHMin: SH=SHq+SHm*dHMin #tiene costanti le altezze per shell con D < dHMin
    
    #mTh=float(minThick)/D
    s.append ((9*"%s\t"+"%s")%("shell#","D","alfa","SH","thick","aColl","ff","fov","ShellWeight(kg)","TotWeight"))    
    while (D>=Dmin):
        #trova diametri con procedura ricorsiva
        alfa=atan2(D,2*FL)/4
        ff=max(fov,alfa/2)
        if (alfa/2)<fov: print n    #debug
        vsh=volume(D,alfa,t,SH)*density/1000000
        v=v+vsh
        #print n,D
        aColl=((D+2*SH*tan(alfa))**2-D**2)*pi/4
        s.append ((9*"%s\t"+"%s")%(n,D,alfa,SH,t,aColl,ff,fov,vsh,v))
        if coating:
            aeff=aColl*(array(coating.reflex(alfa,ener))**2)           
            totAeff=totAeff+aeff
        #passa alla prossima shell
        #t=mTh*D
        #t=float(fixedThick)
        #deltaD=2*SH*tan(alfa)
        D=nextInner(D,fov*60*180/math.pi,(SHm,SHq))        
        SH=SHq+SHm*D
        if dHMin:
            if D<dHMin: SH=SHq+SHm*dHMin
        D=inverso(FL,D-2*t,SH)
        n=n+1
    return "\n".join(s),totAeff

def saw(fileIn,fileOut=None):
    '''crea un file per plottare il profilo delle shell con gnuplot.
    usando with vectors nohead. deve contenere su 4 colonne per ogni shell
    D,0 , deltaD,H'''
    if fileOut==None: fileOut="_saw".join(os.path.splitext(fileIn))
    l=open(fileIn,"r").read().split("\n")
    dah=[(float(ll.split()[1]),float(ll.split()[2]),float(ll.split()[3])) for ll in l[1:]]
    g=open(fileOut,"w")
    for d in dah:
        deltaD=2*d[2]*tan(d[1])
        g.write("%s\t%s\t%s\t%s\n"%(d[0],0,deltaD,d[2]))
        #g.write("%s\t%s\n"%(d[0]+deltaD,d[2]))
    g.close()
    

if __name__=='__main__':
    wfxtDiam(76.,367.158,1600.,((76.,150.),(365.,150.)),0.3,30.,8.8)

    
if __name__=="x__main__":
    '''calcola peso (kg) e dimensioni di un telescopio wfxt con le caratteristiche delle shell
    secondo appunti conconi, vale a dire ff = fov se fov>alfa/2, senno' alfa/2 con fov raggio
    del campo di vista. tutte le grandezze in mm. resituisce stringa con tabella.
    il diametro max e' quello esterno alla pupilla d'ingresso, quello min e' interno all'intersezione'''
    outdir="edge"
    nomeFile="d1571_farf27_t2h"
    ener=[i*0.1 for i in range(1,101)]
    cc=Ccoating ("iridio.dat",ener)
    #(Dmin,Dmax,FL,SHmin,SHmax,fixedThick,fovArcMin,density,coating=None,ener=None)
    #                                                   #DESCRIZIONE            NOME
    #versione 2.0 verifica su un po' di configurazioni vecchie:
    #a,aeff=wfxtDiam(350,700,3000,77.14,120,1.,30,3.3,cc) #F=3m                v1 base_farf30
    #a,aeff=wfxtDiam(350,diretto(3000,708.914695536,121.091668145)+2,3000,77.14,121.091668145,1.,30,3.3,cc) #F=3m                v2 base_farf30
    #a,aeff=wfxtDiam(350,700,3000,77.14,120,2.,30,3.3,cc) #F=3m spess=2        v1 base_farf30_t2
    #a,aeff=wfxtDiam(350,700,3000,77.14,120,2.,30,3.3,cc) #F=3m spess=2        v2 base_farf30_t2
    #a,aeff=wfxtDiam(300,700,2750,55.8,94.286,2.,30,3.3,cc) #new F=2.75m d 300-700 v1 d3070_farf27_t2
    
    #a,aeff=wfxtDiam(300,710,2750,((350,700),(60.612,94.286)),2.,30,3.3,cc) #farf F=2.75m d 300-700 v2 d3070_farf27_t2
    #a,aeff=wfxtDiam(300,710,2750,((350,700),(88,88)),2.,30,3.3,cc) #eq F=2.75m d 300-700 v2 d3070_eq27_t2
    #a,aeff=wfxtDiam(300,710,2750,((350,700),(60.612,94.286)),1.,30,3.3,cc) #farf F=2.75m d 300-700 v2 d3070_farf27_t1
    #a,aeff=wfxtDiam(300,710,2750,((350,700),(90,90)),1.,30,3.3,cc) #eq F=2.75m d 300-700 v2 d3070_eq27_t1

    #a,aeff=wfxtDiam(150,710,2750,((350,710),(60.612,94.286)),2.,30,3.3,cc) #farf F=2.75m d 150-710 v2 d1571_farf27_t2
    #a,aeff=wfxtDiam(150,710,2750,((350,710),(60.612,94.286)),1.,30,3.3,cc) #farf F=2.75m d 50-710 v2 d1571_farf27_t1
    a,aeff=wfxtDiam(150,710,2750,((350,710),(60.612,94.286)),2.,30,3.3,cc,dHMin=300) #farf F=2.75m d 150-710 v2 d1571_farf27_t2h
    open(os.path.join(outdir,"seq_%s.txt"%nomeFile),"w").write(a)
    f=open(os.path.join(outdir,"aeff_%s.txt"%nomeFile),"w")
    for e,a in zip(ener,aeff):f.write(str(e)+"\t"+str(a)+"\n")
    f.close()
    saw(os.path.join(outdir,"seq_%s.txt"%nomeFile))
    plottaFile(os.path.join(outdir,"aeff_%s.txt"%nomeFile),title=nomeFile)
    
def standardGeo():
    from CsimpleCreaDiam import shellGeo
    #(self,focal,shellHeight,Dmax,fillingFactor,nshell,wallDensity,thickVal=XMM)
    #tel=shellGeo(3.5,14,700,0.5,37,3.3,(0,1.,1.)) #geo 35
    tel=shellGeo(3.5,14,700,0.5,37,3.3,(0,1.,1.)) #geo30
    ener=[float(i)/10 for i in range(1,100)]
    o=open(r"edge\Geo_30.txt","w")
    for e,a in zip (ener,tel.aeff(ener)):
        o.write("%s\t%s\n"%(e,a))
    o.close()