
'''todo:
    importando il modulo dice:
   C:\work\pythonAnalyzer\geometria\geoShell.py:1: DeprecationWarning:
   Non-ASCII character '\xb0' in file I:\vince\prog\offaxis\OA\simX2006\generaPlot.py
   on line 211, but no encoding declared;
   see http://www.python.org/peps/pep-0263.html for details 
   forse:
   possibilita' di aggiungere numero iniziale ai nomi file, in modo da averli ordinati
   routine per settare il titolo nei grafici, i valori di focale ecc, nella tabella,
   ed eventualmente il nome dei file di output
istruzioni:
   partendo da elenco di cartelle di risultati (ottenute con i programmi fortran offAxis)
   nel di file di testo in fl,
   per ogni cartella elencata valuta peso, campo di vista interpolato all'energia entarget
   area e angoli usati per il calcolo.
   genera file 'plotta.plt' per il plot con gnuplot di curve di area efficace e dati valutati.
   scrive gli stessi dati in un file 'tabella.txt'
per usarlo:
   impostare enTarget
   impostare in fl l'elenco dei file di testo contenenti gli elenchi di directory
   attenzione alla routine extractTitle che estrae i valori da usare nei titoli dei grafici
   e nelle tabelle dal nome delle directory (nettamente migliorabile, es. potrebbe prenderli
   dal file dei parametri)
note:
   i valori di area efficace e peso indicato nei grafici sono rispettivamente, ridotti del 10%
   (per lo spider) e aumentato del 30% (per la struttura). I valori delle tabelle non comprendono
   correzioni.
   '''

import os
import string
def evaluateFOV(folder,entarget,ang):
    '''valuta il campo di vista all'energia energy per il file delle aree
    offaxis contenuto in folder'''
    
    def createDic(ener,aree,entarget,ang):
        '''data una lista di liste(una per ogni angolo) dell'aree, con le aree in funzione dell'energia
        restituisce per ogni angolo un dizionario con ang l'area interpolata tra l'iesima e l'i+1 esima
        '''
        def findEnTarget(ener,target):
            '''dato un vettore float delle energie ener, restituisce indice di
            quella immediatamente inferiore a target'''
            l= [i for i in ener if i < target]
            ll=len(l)
            return ll-1

        ind= findEnTarget(ener,entarget)
        d={}
        for i in range(len(aree)):
            d[ang[i]]=aree[i][ind]+((aree[i][ind+1]-aree[i][ind])/(ener[i+1]-ener[i]))*(entarget-ener[ind])
        return d
    def interpolateFOV(y,fractionArea=0.5):
        '''dal dizionario {x:y} con x angoli offaxis e y aree normalizzate
        restituisce l'angolo interpolato con area=fractionArea '''
        #print oaDic
        k=y.keys()
        k.sort()
        maxAng=-1
        for x in k:
            if y[x]>fractionArea:
                minAng=(x,y[x])
            elif y[x]<fractionArea:
                maxAng=(x,y[x])
                break
            elif y[x]==fractionArea:
                minAng=maxAng=(x,y[x])
            if maxAng ==-1:        #se > dell'angolo massimo fissa valori per estrapolaizone
                minAng=(k[-2],y[k[-2]])
                maxAng=(k[-1],y[k[-1]])
        if maxAng[0]==minAng[0]: fov = maxAng[0]
        else:fov=(((maxAng[0]-minAng[0])/(maxAng[1]-minAng[1]))*(fractionArea-maxAng[1]))+maxAng[0]
        return fov,fractionArea
    

    areeFile=os.path.join(folder,"aree.txt")
    ener=loadEn(areeFile)    
    #entarget contiene il valore di energia immediatamente inferiore a energy
    #ora serve un pezzo che prende in dizionario angolo:area interpolata
    aree=generalRoutines.loadCol(areeFile,1)
    offAxisAr={}
    offAxisAr=createDic(ener,aree,entarget,ang)
    onAx=offAxisAr[0]
    norm=[(u,v/onAx) for u,v in offAxisAr.items()]
    for i in norm:
        if i[1] >1.0:
            print "casino in aree per E=", entarget
            print "per ",areeFile,":"
            print norm,i
            return -1,-1
    norm=dict(norm)
    #print norm
    FOV,intAr = interpolateFOV(norm)
    return FOV,onAx,intAr*onAx

def extractTitle(st):
    '''dalla stringa con il nome della directory estrae le informazioni
    sulla configurazione e restituisce una lista con
    coating, filling factor, focale, diametro massimo'''
    st=str(st)
    coating={"ir":"Ir monolayer","ml":"W/Si multilayer"}
    ff={"ff000":"0.0","ff00":"0.0","ff007":"0.07","ff015":"0.15"}
    diam={"max60":"60","max70":"70"}
    focale={"20m":"20","22m":"22.5","25m":"25","27m":"27.5","30m":"30"}
    l=[string.lower(w) for w in st.split("_")]
    #print l
    k=[coating[l[0]],ff[l[1]],focale[l[2]],diam[l[3]]]
    #print k
    return k

def offAxisInfo(folder,entarget):
    '''passando un folder e l'energia entarget restituisce
    peso,lista degli angoli, fov interpolato e area corrispondente'''  
    #print folder
    ang=retrieveAngleSteps(folder)
    Weight=Wrecover(folder)
    fov,fovAr=evaluateFOV(folder,entarget,ang)
    return Weight,ang,fov,fovAr

guarda che c è una versione più aggiornata: generaPlot2
if __name__=="__main__":
    entarget=30.
    plt=open("plottaAll.plt","w")
    fl=['tutti.txt']  
    tabella=[]

    for filelista in fl:
        listaDir=open(filelista,"r")
        a=listaDir.read().split()
        for config in a:
            print config
            Weight,ang,FOV,fovAr=offAxisInfo(config ,entarget)
            tabella.append([config ,Weight,FOV,fovAr,ang])
            plt.write("set title ' %s, ang. shell sep.=%s, FL= %s m, Dmax=%s cm' \n"  %tuple(extractTitle(config)))
            plt.write("set term win\n")
            plt.write("set out\n")
            plt.write("set key box\n")
            plt.write("set grid\n")
            plt.write("unset label\n")
            plt.write("set xlabel 'Energy (keV)'\n")
            plt.write("set ylabel 'Effective Area (cm^2) x 90%'\n")
            plt.write('set label "WEIGHT (+30%%): %4.2f Kg\\nFOV@%4.2f keV: %4.2f arcmin" at graph 0.66,0.58\n' %(Weight*1.3,entarget,FOV*2))
            plt.write("plot '" + config + "\\aree.txt' u 1:($2*0.9) index 0 w l title '"+str(ang[0])+" arcmin',\\")
            for i in range(1,len(ang)):
                plt.write("\n'" + config  + "\\aree.txt' u 1:($2*0.9) index "+str(i)+" w l title '"+str(ang[i])+" arcmin',\\")
            plt.write("\n'" + config  + "\\aree.txt' u 1:($2*0.5*0.9) index 0 w l lw 2 title '50% on Axis Area'\n")
            plt.write("\n")    
            plt.write("pause -1\n")
            plt.write("set term png\n")
            plt.write("set output '"+config +".png'\n")
            plt.write("replot\n")
            plt.write("set term win\n")
            plt.write("set out\n")
            plt.write("\n")
        listaDir.close()

    plt.close()

    tab=open("tabella.txt","w")
    tab.write("Coating\tangular shell sep (°) \tFocal lenght(m)\tMax diameter(cm)\t")
    tab.write("Weight(Kg w/o struct)\tFOV(arcmin radius)\tArea(cm^2 @FOVangle)\tSampling Angles(arcmin)\n")
    for el in tabella:
        for e in extractTitle(el[0]): tab.write (e+"\t")
        for dato in el[1:-1]: tab.write(str(dato)+"\t")
        tab.write(str(el[-1])+"\n")
    tab.close()
