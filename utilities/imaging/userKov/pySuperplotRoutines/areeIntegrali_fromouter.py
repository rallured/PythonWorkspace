import os
import string
from userKov.pyGeneralRoutines import generalRoutines
'''dalla piu' esterna alla piu' interna, come nel file'''


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
    dal piu' interno al piu' esterno'''
    ff=open(os.path.join(baseDirW,"pesi.txt"),"r")
    pes=map (float,ff.readlines()[4:])[::-1] #inverte
    ff.close()
    r=[reduce(lambda x, y: x+y, pes[0:i]) for i in range(1,len(pes)+1)]
    '''tt=open("test.txt","w")
    for o in r: tt.write("%s\n"%o)
    tt.close()'''
    return r
    
def recuperaArea(dir,enerTarget=30.0):
    '''recupera l'area efficace dalla directory dir, per l'i-esima shell
    (dalla piu' interna alla piu' esterna)
    all'energia ener, restituisce un vettore conle aree efficaci integrate'''
    nshell=100
    aeffInt=[]
    aT=0
    for i in range(1,nshell+1): #,-1
        areaFile=open(os.path.join(baseDir,'%03i'%(i),'aeff%04i.dat'%(i)),"r")
        ener=[]
        area=[]
        for l in areaFile.readlines()[:-8]:
            ener.append(float(l.split()[0]))
            area.append(float(l.split()[1]))
        areaFile.close()
        aT=aT+interpola(ener,area,[enerTarget])[0]
        aeffInt.append(aT)
    return aeffInt

def recuperaAcoll(dir,dummy):
    '''recupera l'area di raccolta dal file diamplot.txt dalla directory dir, per l'i-esima shell
    (dalla piu' interna alla piu' esterna)
    all'energia ener, restituisce un vettore con le aree di raccolta integrate'''
    nshell=100
    aeffInt=[]
    aT=0
    areaFile=open(os.path.join(baseDir,'diamplot.txt'),"r")
    pes=[a.split()[1] for a in areaFile.readlines()[1:]]
    areaFile.close()
    pes=map (float,pes)        #[::-1]
    r=[reduce(lambda x, y: x+y, pes[0:i]) for i in range(1,len(pes)+1)]
    return r

def recuperaDiam(dir):
    '''recupera diametri dal file diamplot.txt dalla directory dir, per l'i-esima shell
    (dalla piu' interna alla piu' esterna)'''
    nshell=100
    aeffInt=[]
    aT=0
    areaFile=open(os.path.join(baseDir,'diamplot.txt'),"r")
    pes=[a.split()[0] for a in areaFile.readlines()[1:]]
    areaFile.close()
    pes=map (float,pes)        #[::-1]
    return pes

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

    diametri=recuperaDiam(dir)    
    for quantita in listaQuant:
        #calcola vettore quantita' integrate per shell
        qtot=quantita[1](dir,quantita[2])

        #crea file di output
        outName="%s_%s.txt"%(quantita[0],f)
        outData=open(os.path.join(resdir,outName),"w")
        outData.write("nShell\tDiam\t%s\t%s x %s\n"%(quantita[0],quantita[0],quantita[4]))
        factor=quantita[4]
        for i,q in enumerate(qtot):
            outData.write("%s\t%s\t%s\t%s\n"%(i+1,diametri[i],q*factor,q))
        outData.close()

        #crea file di plot
        conf="%s, F=%s, ff=%s, Dmax=%s" %tuple(extractTitle(f))
        plt.write("set ylabel '%s'\n" %quantita[3][0])
        plt.write("set xlabel 'Shell Number from inner'\n" )
        plt.write("set title '%s for %s configuration'\n" %('Integrated '+"".join(quantita[3]),conf))
        plt.write("plot '%s' u 1:2 title '%s x %s' w lp\n" %(outName,quantita[3][0],quantita[4]))
        plt.write("pause-1\n\n")

    '''plotta insieme i valori degli elementi di listaQuant indicati da sel,
    usando come etichetta '''    
    sel=[[0,1,2,3,4,5,"[0:1200]"],[6,"[0:50]"]]
    #pltTot.write("unset key\n")
    pltTot.write("set title '%s'\n"%conf)
    pltTot.write("set key top left box\n")
    pltTot.write("set xlabel 'Shell Number from inner'\n")
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

if __name__=="x__main__":
    #i parametri sono impostati ora,
    #si potrebbero far recuperare dalle impostazioni
    resdir="HXMT_8_integrali_test"
    nshell=100
    wdName="HXMT_8"
    offAxisDir="G:\\vince\\prog\\f_2005\\ufficiali\\prog_superplotGen1\\HXMT_8\\hxmtStart\\"
    dirfomName="100Sh_dmax40"
    f="testNshell_plt"
    baseDir=os.path.join(os.path.pardir,wdName,f,dirfomName)
    baseDirW=os.path.join(offAxisDir,"100Sh_dmax40")
    listaQuant=[['area02',recuperaArea,2.,('Effective Area (cm^2) ','@ 2 keV'),1.],
            ['area03',recuperaArea,3.,('Effective Area (cm^2) ','@ 3 keV'),1.],
            ['area04',recuperaArea,4.,('Effective Area (cm^2) ','@ 4 keV'),1.],
            ['area06',recuperaArea,6.,('Effective Area (cm^2) ','@ 6 keV'),1.],
            ['area08',recuperaArea,8.,('Effective Area (cm^2) ','@ 8 keV'),1.],
            ['Acoll',recuperaAcoll,0,('Collecting Area (cm^2) ',''),1.],
            ['weight',recuperaPeso,"",('Weight (kg)',''),1.]]
    #elenco di liste contenenti nome di quantita e routine per recuperare i valori
    #nel formato[filePrefix,routine,parametro,(labely,titleSuffix),structFactor]
    HXMT()
    
if __name__=="__main__":
    #i parametri sono impostati ora,
    #si potrebbero far recuperare dalle impostazioni
    resdir="HXMT_8_integrali_test"
    nshell=100
    wdName="HXMT_8"
    offAxisDir="G:\\vince\\prog\\f_2005\\ufficiali\\prog_superplotGen1\\HXMT_8\\hxmtStart\\"
    dirfomName="100Sh_dmax40"
    f="testNshell_plt"
    baseDir=os.path.join(os.path.pardir,wdName,f,dirfomName)
    baseDirW=os.path.join(offAxisDir,"100Sh_dmax40")
    listaQuant=[['area02',recuperaArea,2.,('Effective Area (cm^2) ','@ 2 keV'),1.],
            ['area03',recuperaArea,3.,('Effective Area (cm^2) ','@ 3 keV'),1.],
            ['area04',recuperaArea,4.,('Effective Area (cm^2) ','@ 4 keV'),1.],
            ['area06',recuperaArea,6.,('Effective Area (cm^2) ','@ 6 keV'),1.],
            ['area08',recuperaArea,8.,('Effective Area (cm^2) ','@ 8 keV'),1.],
            ['Acoll',recuperaAcoll,0,('Collecting Area (cm^2) ',''),1.],
            ['weight',recuperaPeso,"",('Weight (kg)',''),1.]]
    #elenco di liste contenenti nome di quantita e routine per recuperare i valori
    #nel formato[filePrefix,routine,parametro,(labely,titleSuffix),structFactor]
    HXMT()