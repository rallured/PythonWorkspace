import os
import string
from userKov.pyGeneralRoutines import generalRoutines
from numpy import interp
'''dalla piu' interna alla piu' esterna, cioe' al contrario del file'''


def extractTitle(st): #da plotRoutines
    '''dalla stringa con il nome della directory estrae le informazioni
    sulla configurazione e restituisce una lista con
    coating, filling factor, focale, diametro massimo'''
    st=str(st)
    coating={"ir":"Ir monolayer","ml":"W/Si multilayer","m2":"Pt/C multilayer"}
    ff={"ff000":"0.0","ff00":"0.0","ff007":"0.07","ff015":"0.15"}
    diam={"max60":"60","max70":"70"}
    focale={"20m":"20","22m":"22.5","25m":"25","27m":"27.5","30m":"30"}
    l=[string.lower(w) for w in st.split("_")]
    #print l
    try:
        k=[coating[l[0]],ff[l[1]],focale[l[2]],diam[l[3]]]
    #print k
    except:
        return [st,"","",""]
    return k

def recuperaPeso(dir,enerTarget=""):
    '''recupera il peso delle shell dalla directory dei risultati offAxis,
    sono in ordine dal piu' esterno al piu' interno'''
    ff=open(os.path.join(baseDirW,"pesi.txt"),"r")
    pes=map (float,ff.readlines()[4:])
    ff.close()
    r=[reduce(lambda x, y: x+y, pes[0:i]) for i in range(1,len(pes)+1)]
    '''tt=open("test.txt","w")
    for o in r: tt.write("%s\n"%o)
    tt.close()'''
    return r
    
def recuperaArea(dir,enerTarget=30.0):
    '''recupera l'area efficace dalla directory dir, per l'i-esima shell
    (dalla piu' esterna alla piu' interna, nel file sono al contrario)
    all'energia ener, restituisce un vettore conle aree efficaci integrate'''
    aeffInt=[]
    aT=0
    for i in range(nshell,0,-1):
        areaFile=open(os.path.join(baseDir,'%03i'%(i),'aeff%04i.dat'%(i)),"r")
        ener=[]
        area=[]
        for l in areaFile.readlines()[:-8]:
            ener.append(float(l.split()[0]))
            area.append(float(l.split()[1]))
        areaFile.close()
        aT=aT+interp([enerTarget],ener,area)[0]
        aeffInt.append(aT)
    return aeffInt

def recuperaDiam(folder):
    '''recupera l'area di raccolta dal file diamplot.txt dalla directory dir, per l'i-esima shell
    (dalla piu' esterna alla piu' interna, nel file sono al contrario)
    all'energia ener, restituisce un vettore con le aree di raccolta integrate'''
    
    areaFile=open(os.path.join(folder,'diamplot.txt'),"r")
    lines=areaFile.readlines()[1:]
    areaFile.close()
    diam=[a.split()[0] for a in lines]
    diam=map (float,diam)[::-1]
    return diam

def recuperaAcoll(dir,dummy):
    '''recupera l'area di raccolta dal file diamplot.txt dalla directory dir, per l'i-esima shell
    (dalla piu' esterna alla piu' interna, nel file sono al contrario)
    all'energia ener, restituisce un vettore con le aree di raccolta integrate'''

    areaFile=open(os.path.join(baseDir,'diamplot.txt'),"r")
    lines=areaFile.readlines()[1:]
    acoll=[a.split()[1] for a in lines]
    areaFile.close()
    acoll=map (float,acoll)[::-1]
    r=[reduce(lambda x, y: x+y, acoll[0:i]) for i in range(1,len(acoll)+1)]
    return r

def files(lista):
    ff=open("interessanti.txt","r")
    l=ff.read().split()
    ff.close()
    return l
    #return ['m2_ff015_20m_max70'] #per test

def replot (pltfile,imgFile,term):
    '''scrive sul file pltfile l'struzione per un replot sul terminale term
    per generare con gnuplot il file grafico di nome imgFile'''
    pltfile.write("set terminal %s\n" %term)
    pltfile.write("set output '%s'\n" %imgFile)
    pltfile.write("replot\n")
    pltfile.write("set term win\n")
    pltfile.write("set out\n\n" )
    

def HXMT():
    if os.path.exists(resdir):
        print "la directory %s esiste gia'!\n" %resdir
        yn=raw_input('sovrascrivo? (y/n)')
        if yn<>"y" : exit
    else: os.mkdir(resdir)
    
    plt=open(os.path.join(resdir,"separati.plt"),'w')
    pltTot=open(os.path.join(resdir,"riepilogo.plt"),"w")
    diam=recuperaDiam(baseDir)
    
    for quantita in listaQuant:
        #calcola vettore quantita' integrate per shell
        qtot=quantita[1](dir,quantita[2])

        #crea file di output e scrive intestazione
        outName="%s_%s.txt"%(quantita[0],f)
        outData=open(os.path.join(resdir,outName),"w")
        outData.write("nShell\tDiam\t%s\n"%quantita[0])
        factor=quantita[4]
        for i,(d,q) in enumerate(zip(diam,qtot)):
            outData.write("%s\t%s\t%s\t%s\n"%(i+1,d,q*factor,q))
        outData.close()

        #crea file di plot
        conf="%s, F=%s, ff=%s, Dmax=%s" %tuple(extractTitle(f))
        plt.write("set ylabel '%s'\n" %quantita[3][0])
        plt.write("set xlabel 'Shell Number from outer'\n" )
        plt.write("set title '%s for %s configuration'\n" %('Integrated '+"".join(quantita[3]),conf))
        plt.write("plot '%s' u 1:3 title '%s x %s' w lp\n" %(outName,quantita[3][0],quantita[4]))
        plt.write("pause-1\n\n")

    '''plotta insieme i valori degli elementi di listaQuant indicati da sel,
    usando come etichetta '''    
    #pltTot.write("unset key\n")
    pltTot.write("set title '%s'\n"%conf)
    pltTot.write("set key top left box\n")
    pltTot.write("set xlabel 'Shell Number from outer'\n")
    for i,trackGroup in enumerate(sel):
        pltTot.write("set yrange %s \n"%trackGroup[-1])
        pltTot.write("set ylabel '%s'\n"%(listaQuant[trackGroup[0]][3][0]))
        pltTot.write("plot ")
        for gg in trackGroup[0:-2]:
            quantita=listaQuant[gg]
            outName="%s_%s.txt"%(quantita[0],f)
            pltTot.write("'%s' u 1:3 title '%s x %s' w lp,\\\n"%(outName,"".join(quantita[3]),listaQuant[trackGroup[0]][4]))
        quantita=listaQuant[trackGroup[-2]]
        outName="%s_%s.txt"%(quantita[0],f)
        pltTot.write("'%s' u 1:3 title '%s x %s' w lp\n"%(outName,"".join(quantita[3]),listaQuant[trackGroup[0]][4]))
        pltTot.write("pause-1\n\n")
        replot(pltTot,"%03i_"%i+os.path.splitext(outName)[0]+".ps","postscript color")
        replot(pltTot,"%03i_"%i+os.path.splitext(outName)[0]+".png","png")
    pltTot.close()
    plt.close()

if __name__=="__main__":

    '''    
    #i parametri sono impostati ora,
    #si potrebbero far recuperare dalle impostazioni
    resdir="HXMT_8_integrali_test"
    nshell=100
    wdName="HXMT_8"
    offAxisDir=r'H:\vince\risultati\HXMT\dati\HXMT_8\hxmtStart_focals\\'
    dirfomName="100Sh_dmax40_f245"
    f="testNshell_plt"
    baseDir=os.path.join(os.path.pardir,wdName,f,dirfomName)
    baseDirW=os.path.join(offAxisDir,"100Sh_dmax40")
    listaQuant=[['area02',recuperaArea,2.,('Effective Area (cm^2) ','@ 2 keV'),1.],
            ['area04',recuperaArea,4.,('Effective Area (cm^2) ','@ 4 keV'),1.],
            ['area06',recuperaArea,6.,('Effective Area (cm^2) ','@ 6 keV'),1.],
            ['area08',recuperaArea,8.,('Effective Area (cm^2) ','@ 8 keV'),1.], 
            ['Acoll',recuperaAcoll,0,('Colleting Area (cm^2) ',''),1.], 
            ['weight',recuperaPeso,"",('Weight (kg)',''),1.]]
    #elenco di liste contenenti nome di quantita e routine per recuperare i valori
    #nel formato[filePrefix,routine,parametro,(labely,titleSuffix),structFactor]
    HXMT()
    '''

    #i parametri sono impostati ora,
    #si potrebbero far recuperare dalle impostazioni
    resdir="HEXITsat_integrali"
    nshell=120
    wdName=r"F:\Nuova\src_plt\hexitSat2009\fullsizeD450_F10_1plt"
    offAxisDir=r'F:\Nuova\src_plt\hexitSat2009\fullsizeD450_F10_1plt'  #per recuperare pesi
    dirfomName="upperlimit"
    f="F10D350ff010_thsx"
    baseDir=os.path.join(wdName,dirfomName)
    baseDirW=os.path.join(offAxisDir)
    sel=[[0,1,2,3,4,5,6,"[0:1000]"],[7,"[0:100]"]]
    listaQuant=[['area10',recuperaArea,10.,('Effective Area (cm^2) ','@ 10 keV'),1.],
            ['area20',recuperaArea,10.,('Effective Area (cm^2) ','@ 10 keV'),1.],
            ['area30',recuperaArea,30.,('Effective Area (cm^2) ','@ 30 keV'),1.],
            ['area50',recuperaArea,50.,('Effective Area (cm^2) ','@ 50 keV'),1.],
            ['area60',recuperaArea,60.,('Effective Area (cm^2) ','@ 60 keV'),1.],
            ['area70',recuperaArea,70.,('Effective Area (cm^2) ','@ 70 keV'),1.],
            ['Acoll',recuperaAcoll,0,('Colleting Area (cm^2) ',''),1.], 
            ['Weight',recuperaPeso,"",('Weight (kg)',''),1.]]
    #elenco di liste contenenti nome di quantita e routine per recuperare i valori
    #nel formato[filePrefix,routine,parametro,(labely,titleSuffix),structFactor]
    HXMT()