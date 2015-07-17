import namelist_class
import os
import math

class shellGeo:
    '''classe che regola il dimensionamento delle shell, a partire
    dai parametri. '''
    XMM=(0.0016,-0.02,0)
    SX=(0.0006,-0.04,0)
    HXMT=(0.0006,0,0.2)
    HXMT_light=(0.0006,0,0.1)

    def diam (self):
        #tutte le misure in cm,
        #risultati convertiti in mm all'uscita
        ind= range(self.nshell)
        m = self.thickVal[0]
        q = self.thickVal[1]
        ds=[]
        spes=[]
        d=self.Dmax
        if d==60:ee=0.3
        if d==70:ee=0.36
        for i in ind:
            s=ee*d/self.Dmax
            if s<self.thickVal[2]:s=self.thickVal[2]
            spes.append(s) #la lista degli spessori e' in mm
            #print (d*10*m + q),spes[-1]
            d=d-2*(spes[-1]/10)
            ang= math.atan2(d,(2*self.FL*100))/4
            dmed=d-2*self.SH*math.tan(ang)
            ds.append(dmed*10)
            d=dmed-2*self.SH*math.sin(math.radians(self.ff))
        return ds, spes
    
    def diam2 (self):
        #tutte le misure in mm,imita fortran
        self.angle=[]
        ind= range(self.nshell)
        m = self.thickVal[0]
        q = self.thickVal[1]
        ds=[]
        spes=[]
        r=self.Dmax/2
        SHmm=self.SH*10
        FLmm=self.FL*1000
        alkw=SHmm*math.sin(math.radians(self.ff))
        for i in ind:
            s=2*r*m + q
            if s<self.thickVal[2]:s=self.thickVal[2]
            spes.append(s)
            r=r-(spes[-1])
            th= math.atan2(r,FLmm)/4
            self.angle.append(th)
            a=2*math.tan(th)*SHmm
            r0=(math.sqrt(a**2+4*r**2)-a)/2
            while True:
                th= math.atan2(r0,FLmm)/4
                a=2*math.tan(th)*SHmm
                r1=(math.sqrt(a**2+4*r**2)-a)/2
                #print i,abs(r1-r0)
                if (abs(r1-r0)<=0.0001):break
                r0=r1
            rmed=r1
            ds.append(rmed*2)
            r=rmed-alkw
        return ds, spes

    def weights(self,diam,spes):
        ang=[math.atan(d/(2000*self.FL))/4 for d in diam]
        W=[]
        for i in range(len(diam)):
            Dmed=diam[i]
            Dmax=Dmed+2*math.tan(ang[i])*self.SH*10
            Dmin=Dmed-2*math.tan(3*ang[i])*self.SH*10
            p=((Dmax+Dmin)/2+Dmed)*math.pi*spes[i]*self.SH*self.wd/100000
            #print Dmax,Dmed,Dmin,spes[i]
            W.append(p)
        return W

    def __init__(self,focal,shellHeight,Dmax,fillingFactor,nshell,wallDensity,thickVal=XMM):
        self.FL=focal
        self.SH=shellHeight
        self.Dmax=Dmax
        self.ff=fillingFactor
        self.nshell=nshell
        self.wd=wallDensity
        self.thickVal=thickVal

    def fill (self,Dmin):
        '''returns the number of shell to fill the diamters larger than dmin, with
        the current configuration parameters by self'''
        store =self.nshell
        self.nshell=50
        n=self.nshell
        while n==self.nshell:
            self.nshell=self.nshell+1
            n=[j>=Dmin for j in self.diam2()[0]].count(True)            
        res=n
        self.nshell=store
        return res


       
class CSimpleCreaDiam(namelist_class.Namelist,shellGeo):
    '''classe per creazione di diametri con parametri letti da namelist,si inizializza con il nome del
    file delle impostazioni. manca qualche perfezionamento
    per gestire tutti i casi della namelist, come i vari flag o il diametro max non fornito direttamente'''
    ##attenzione - le variabili di shellGeo e Namelist non sono sincronizzate, quando cambio i valori di una
    ##non si rispecchiano sull'altra.
    
    def __init__(self,nl):
        namelist_class.Namelist.__init__(self,nl)
        self.dirorigin=os.path.join(self["WORKDIR"].strip("\"' "),self["DIRORIGIN"].strip("\"' "))
        self.dirrisult=os.path.join(self["WORKDIR"].strip("\"' "),self["DIRRISULT"].strip("\"' "))
        shellGeo.__init__(self,self["F_LENGTHdaImp_m"],
                          self["F_HEIGHTdaImp_cm"],self["massimo"],self["FieldOfViewDeg"],
                          self["nshDaImp"],self["WallDensity"],(self["ThickM"],self["ThickQ"],self["minThickmm"]))

class CSimpleCreaDiamOA(namelist_class.Namelist,shellGeo):
    '''come la precedente, con la differenza che nelle namelist non ci sono ThickM, ThikQ, minThickmm
    , ma linearThick e thickExtSh, come nelle vecchie versioni dei programmi OA'''
    
    def __init__(self,nl):
        namelist_class.Namelist.__init__(self,nl)
        self.dirorigin=os.path.join(self["WORKDIR"].strip("\"' "),self["DIRORIGIN"].strip("\"' "))
        self.dirrisult=os.path.join(self["WORKDIR"].strip("\"' "),self["DIRRISULT"].strip("\"' "))
        shellGeo.__init__(self,self["F_LENGTHdaImp_m"],
                          self["F_HEIGHTdaImp_cm"],self["massimo"],self["FieldOfViewDeg"],
                          self["nshDaImp"],self["WallDensity"],(self["ThickExtSh"]/self["massimo"],0,0))

        
        
        
        
if __name__=="__main__":
    nl2=CSimpleCreaDiam("imp_OffAxis.txt")

    focal=[i*2.5 for i in range(8,12)]
    ff=[0.00,0.07,0.15]
    D=[60,70]
    comblist=[[f1,f2,f3] for f1 in focal for f2 in ff for f3 in D]

    #comblist=[[20.0,0.0,60]]
    p=open("confrDiam.plt","w")
    #p2=open("confrDiam2.plt","w")
    #p2.write("plot ")
    for comb in comblist:
        ndf="F"+str(comb[0]).split(".")[0]+"_ff"+str(comb[1]).split(".")[1].zfill(3)+"_D"+str(comb[2])+".txt"
        nl2.FL,nl2.ff,nl2.Dmax=comb
        di,sp=nl2.diam2()
        f= open(ndf,"w")
        for d in di[-1:0:-1]:
            f.write(str(d)+"\n")
        print ndf,nl2.FL,nl2.SH,nl2.Dmax,nl2.ff,nl2.nshell,nl2.wd,(nl2.thickVal[0],nl2.thickVal[1])
        print "\n"
        f.close()
        #per plottare confronto con vecchi diametri
        p.write("plot '"+ndf+"' u 0:1 w lp")
        fcnf1="mlff"+str(comb[1]).split(".")[1].zfill(3)+str(comb[0]).split(".")[0]+"mmax"+str(comb[2])+".txt"
        fcnf2="m2ff"+str(comb[1]).split(".")[1].zfill(3)+str(comb[0]).split(".")[0]+"mmax"+str(comb[2])+".txt"
        if os.path.exists(os.path.join("af_files",fcnf1)):
           p.write(",\\\n'"+os.path.join("af_files",fcnf1)+"' u 0:1 w lp")
        if os.path.exists(os.path.join("af_files",fcnf2)):
           p.write(",\\\n'"+os.path.join("af_files",fcnf2)+"' u 0:1 w lp")
        '''if os.path.exists(os.path.join("frlike",ndf)):
           p.write(",\\\n'"+os.path.join("frlike",ndf)+"' u 0:1 w lp")
        fcnf3="ml_ff"+str(comb[1]).split(".")[1].zfill(3)+"_"+str(comb[0]).split(".")[0]+"m_max"+str(comb[2])+".txt"
        if os.path.exists(os.path.join("extractedDiam",fcnf3)):
           p.write(",\\\n'"+os.path.join("extractedDiam",fcnf3)+"' u 0:1 w lp")
        fcnf4="ir_ff"+str(comb[1]).split(".")[1].zfill(3)+"_"+str(comb[0]).split(".")[0]+"m_max"+str(comb[2])+".txt"
        if os.path.exists(os.path.join("extractedDiam",fcnf4)):
           p.write(",\\\n'"+os.path.join("extractedDiam",fcnf4)+"' u 0:1 w lp")'''
        p.write("\n")
        p.write("pause -1\n\n")

        p.write("set term png\n")
        p.write("set out '"+ndf+".png'\n")
        p.write("replot\n")
        p.write("set term win\n set out\n\n")

        #p2.write("'"+ndf+"' u 0:1 w lp,\\\n")    

    p.close()    
    #p2.close()        
    


if __name__=="x__main__":
    '''questa fa il calcolo di quante shell ci vogliono per coprire un
    range di diametri per tutte le configurazioni lette da una serie di dirfom'''

    path=r"F:\vince\risultati\HXMT\dati\HXMT_11\hxmtStart_focals"
    file="imp_offAxis.txt"

    f=open(os.path.join(path,"dirlist.dat"),"r").read().split()
    Dmin=90
    for i in f:
        a=CSimpleCreaDiam(os.path.join(path,i,file))
        a.nshell=60
        print i, [j>=Dmin for j in a.diam2()[0]].count(True)


def printPeso(spesLaw=CSimpleCreaDiam.HXMT,conf=None,NMaxSh=50):
    '''questa invece lo calcola leggendo una configurazione e variandone una serie di
    parametri passati in conf (lista di liste) con legge di spessori e il massimo numero
    di shell. restituisce stringa con riassunto'''

    if conf==None: conf=[[2.45,30,300,30],[2.65,30,300,30],[2.85,30,300,30],
          [2.45,30,330,30],[2.65,30,330,30],[2.85,30,330,30],
          [2.45,30,350,30],[2.65,30,350,30],[2.85,30,350,30]]    
    
    def setShell(a,v):
            F,H,Dmax,nsh=v
            a.FL=F
            a.SH=H
            a.Dmax=Dmax
            a.nshell=nsh

    def findNshell(a,Dmin,Nmax):
        '''calcolo di quante shell ci vogliono per coprire un
        range di diametri'''
        if Nmax!=-1:
            n=Nmax
            ntry=Nmax
            while (n==ntry and ntry<=Nmax):
                ntry=ntry+1
                a.nshell=ntry
                d,s=a.diam2()        
                Dmin=90
                n = [j>=Dmin for j in d].count(True)
            return min(n,Nmax)
        else: return a.nshell
        
    def printDrange(a):
        rv=[]
        d,s=a.diam2()
        w=a.weights(d,s)
        rv.append( "F=%s, H=%s, Dmax=%s, nshell=%s: "%(a.FL,a.SH,a.Dmax,a.nshell))
        rv.append( "diametri da %s a %s"%(d[-1],d[0]))
        rv.append( "spessori da %s a %s"%(s[-1],s[0]))
        rv.append( "Pesi: da %s a %s"%(w[-1],w[0]))
        rv.append( "Totale: %s\n"%reduce(lambda x, y: x+y, w))
        return "\n".join(rv)

    file=r"F:\vince\risultati\HXMT\dati\hxmt_9\hxmtStart_focals\f245\imp_offAxis.txt"
    a=CSimpleCreaDiam(file)

    Dmin=90
    Nmax=NMaxSh
    strv=[]
    for v in conf:
        setShell(a,v)
        a.thickVal=(spesLaw)
        ns=findNshell(a,Dmin,Nmax)
        strv.append( "Dmin=%s, ci vogliono %s shell,"%(Dmin,ns))
        strv.append( "infatti:")
        a.nshell=ns-1
        strv.append(printDrange(a))
    return "\n".join(strv)

if __name__=="x__main__":
    open("HXMT_py_bh.txt","w").write( printPeso(spesLaw=CSimpleCreaDiam.HXMT,conf=None,NMaxSh=30))
    open("HXMT_py_bl.txt","w").write( printPeso(spesLaw=CSimpleCreaDiam.HXMT_light,conf=None,NMaxSh=30))
    open("HXMT_py_eh.txt","w").write( printPeso(spesLaw=CSimpleCreaDiam.HXMT,conf=None,NMaxSh=50))
    open("HXMT_py_el.txt","w").write( printPeso(spesLaw=CSimpleCreaDiam.HXMT_light,conf=None,NMaxSh=50))

if __name__=="x__main__":
          
    file="G:\\vince\\risultati\\ottimizzazioni\\HXMT\\HXMT_8\\hxmtStart_focals\\100Sh_dmax40_f245\\imp_offAxis.txt"
    a=CSimpleCreaDiam(file)

    a.thickVal=a.HXMT_light    
    dl=a.diam2()
    wl=a.weights(dl[0],dl[1])
    a.thickVal=a.HXMT    
    dh=a.diam2()   
    wh=a.weights(dh[0],dh[1])

    #inverte e integra
    dl=dl[0][::-1]
    wl=wl[::-1]
    wh=wh[::-1]
    wl=[reduce(lambda x, y: x+y, wl[0:i]) for i in range(1,len(wl)+1)]
    wh=[reduce(lambda x, y: x+y, wh[0:i]) for i in range(1,len(wh)+1)]
    
    #output    
    ff=open("Wtrend.dat","w")    
    for i in range(len(wl)):
        ff.write("%s\t%s\t%s\t%s\n"%(i+1,dl[i],wl[i],wh[i]))
    ff.close()


def NshellPerD (listaFile):
    '''analizza una serie di configurazioni a partire dai file delle impostazioni contenuti nelle
    cartelle elencate nel file listaFile.
    ne modifica il diametro massimo e calcola il numero di shell per lo stesso diametro minimo.
    genera i nuovi file di impostazioni per il raytracing. serve per raggruppare piu' configurazioni
    in un singolo raytracing per avere plot incrementale'''

    outDir="inc_conf"
    confs=open(listaFile,"r").readlines()
    if not os.path.isdir(outDir): os.mkdir(outDir)
    count=1
    for c in confs:
        confName,Dmin=c.split()[0:2]
        wolter=CSimpleCreaDiamOA(os.path.join(confName,"imp_OffAxis.txt"))
        wolter.Dmax=700
        wolter["massimo"]=700
        wolter["nshDaImp"]=wolter.fill(float(Dmin))
        wolter.write(os.path.join(outDir,confName+"_imp.txt"))
        open(os.path.join(outDir,confName+"_imp.txt"),"a").write(str(wolter["nshDaImp"]))
        print "scritti i risultati per la directory %s,\nanalizzati %s cartelle" %(confName,count)
        print "servono %s shell"%wolter["nshDaImp"]
        count=count+1

if __name__=="x__main__":
    os.chdir(r"C:\work\copiaOA\batch\simx2007\\")
    NshellPerD(r"dminMax.txt")

if __name__=="x__main__":
    l=os.listdir(r"C:\work\copiaOA\simX2006\M2imp_files")
    for ll in l[0:1]:
        a=CSimpleCreaDiamOA(ll)
        a["dirrisult"]=a["dirrisult"].lower().replace("ml","m2")
        a.write(r"C:\work\copiaOA\simX2006\M2imp_files\m2"+"\\%s"%ll)

if __name__=="x__main__":
    os.chdir(r"C:\work\copiaOA\batch\batch4")
    lista=["m2_ff007_195_max65","m2_ff015_185_max65",
           "m2_ff015_190_max65","m2_ff015_190_max70",
           "m2_ff015_200_max65","m2_ff015_200_max70"]

    for name in lista:
        a=CSimpleCreaDiamOA(os.path.join("imp_files",name+".txt"))
        open(os.path.join("imp_files","diam",name+".txt"),"w").write( "\n".join(map(str,a.diam2()[0][::-1])))
        
if __name__=="__main__":
    os.chdir(r"E:\work\WTDf")
    a=CSimpleCreaDiamOA("imp_OffAxis.dat")
    