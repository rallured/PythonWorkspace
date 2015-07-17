from math import *
def findDMed(Dmax, F=20000.,H=300.,mode=0,toll=0.0001):
    '''crea ricorsivamente il diametro all'intersezione par-hyp, partendo dal diametro
    interno della pupilla di entrata. F focale, H lunghezza shell. per confronto sono mostrati
    anche i valori ricavati con il sistema di S.Basso, cioe' approssimando tan e atan con l'argomento.
    mode indica con che metodo effettuare il calcolo per il risultato restituito (per default Conconi)
    se < 0 restituisce una lista [conconi, basso, kov]'''

    r=Dmax/2
    th= atan2(r,F)/4
    a=2*tan(th)*H
    r0=(sqrt(a**2+4*r**2)-a)/2
    r0k=r-H*tan(th)
    rB=r*(4*F)/(H+4*F)
    thB=atan2(rB,F)/4
    for i in range (10):
        th= atan2(r0,F)/4
        thk=atan2(r0k,F)/4
        a=2*tan(th)*H
        r1=(sqrt(a**2+4*r**2)-a)/2
        r1k=r-H*tan(thk)
        if mode<0:print "\n-%s-----------\n\tConconi\tBasso\tkov"%(i)
        if mode<0:print "MAX: %s\t%s\t%s" %(r0+tan(th)*H,rB+tan(thB)*H,r0k+tan(thk)*H)
        if mode<0:print "%s: %s\t%s\t%s \n%s: %s\t%s \n%s: %s\t%s" %("r0",r0,rB,r0k,"r1",r1,r1k,"difference",r1-r0,r1k-r0k)
        if mode<0:print "angoli: %s\t%s\n---------\n"%(th,thB)
        r0=r1
        r0k=r1k
    if mode<0: return [r0*2,rB*2,r0k*2]
    else: return [r0*2,rB*2,r0k*2][mode]

def findNShell (Dmax,Dmin,spes,ff, F=20000.,H=300.,mode=0,toll=0.0001):
    '''restituisce il numero di shell con filling factor ff,
    trascurando spessori shell
    dati diam massimo e minimo e le impostazioni da passare a findDmed'''
    d=Dmax-2*spes*cos(atan(Dmax/(8*F)))
    n=0
    while (d>Dmin):
        d=findDMed(d, F,H,mode,toll)-2*H*tan(ff*180/pi)-2*spes*cos(atan(d/(8*F)))
        print d
        n=n+1
    return n-1
        
    


if __name__=="__main__":
    dStart=700
    #print findDMed(float(dStart))
    print findNShell(300,115,0.2,0,2450,300)