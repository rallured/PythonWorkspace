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

def wfxtDiamMix(Dmin,Dmax,FL,SHpoints,fixedThick,fovArcMin,density,coatIn=None,coatOut=None,i=None,ener=None,dHMin=None):
    '''calcola sequenza di diametri. SHpoints sono due coppie di coordinate in grafico (D,SH). '''
    s=[] #sequenza di stringhe che vengono unite alla fine e restituite come tabella

    fov=radians(float(fovArcMin)/60)
    v=0
    n=1
    totAeff=0
    t=float(fixedThick)    
    SHm,SHq=rPars(SHpoints)
    if dHMin:
        if Dmax<dHMin: SH=SHq+SHm*dHMin #tiene costanti le altezze per shell con D < dHMin
    SH= SHm*Dmax+SHq  #float(SHmin)
    D=inverso(FL,Dmax-2*t,SH)  #float(Dmin)
    j=0
    #mTh=float(minThick)/D
    s.append ((9*"%s\t"+"%s")%("shell#","D","alfa","SH","thick","aColl","ff","fov","ShellWeight(kg)","TotWeight"))    
    while (D>=Dmin):
        j=j+1
        #trova diametri con procedura ricorsiva
        alfa=atan2(D,2*FL)/4
        ff=max(fov,alfa/2)
        #if (alfa/2)<fov: print n    #debug
        vsh=volume(D,alfa,t,SH)*density/1000000
        v=v+vsh
        aColl=((D+2*SH*tan(alfa))**2-D**2)*pi/4
        s.append ((9*"%s\t"+"%s")%(n,D,alfa,SH,t,aColl,ff,fov,vsh,v))
        if j<=i:
            coating=coatOut
        else:
            coating=coatIn
        #print n,coating
        aeff=aColl*(array(coating.reflex(alfa,ener))**2)           
        totAeff=totAeff+aeff
        #passa alla prossima shell
        #t=mTh*D
        #t=float(fixedThick)
        #deltaD=2*SH*tan(alfa)
        if dHMin:
            if D<dHMin: SH=SHq+SHm*dHMin
        D=nextInner(D,fov*60*180/math.pi,(SHm,SHq))        
        SH=SHq+SHm*D
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

if __name__=="__main__":
    '''prova  a combinare shell in nichel con shell in ir
    calcola peso (kg) e dimensioni di un telescopio wfxt con le caratteristiche delle shell
    secondo appunti conconi, vale a dire ff = fov se fov>alfa/2, senno' alfa/2 con fov raggio
    del campo di vista. tutte le grandezze in mm. resituisce stringa con tabella.
    il diametro max e' quello esterno alla pupilla d'ingresso, quello min e' interno all'intersezione'''
   
    outdir="edge\\edgeMix"
    nomeFile="d2175_farf27_t1Nimix"
    ener=[i*0.1 for i in range(1,101)]
    cc=Ccoating ("iridio.dat",ener)
    cn=Ccoating ("Nichel.dat",ener)
    
    plt=open(os.path.join(outdir,"plotMix%s.plt"%nomeFile),"w")

    for n in range(68):  #68
        a,aeff=wfxtDiamMix(210,750,2750,((350,700),(60.612,94.286)),2.,30,3.3,cn,cc,n) #farf F=2.75m d 210-750 v2 d2175_farf27_t2Ni    
        plt.write ("plot 'aeff_%s%s.txt' u 1:2 w l\npause-1\n\n "%(nomeFile,n))

        f=open(os.path.join(outdir,"aeff_%s%s.txt"%(nomeFile,n)),"w")
        for e,a in zip(ener,aeff):f.write(str(e)+"\t"+str(a)+"\n")
        f.close()

        #plottaFile(os.path.join(outdir,"aeff_%s%s.txt"%nomeFile),title=nomeFile)
        
    #open(os.path.join(outdir,"seq_%s.txt"%(nomeFile)),"w").write(a)
    #saw(os.path.join(outdir,"seq_%s.txt"%nomeFile))  
    plt.close() 
        
        
