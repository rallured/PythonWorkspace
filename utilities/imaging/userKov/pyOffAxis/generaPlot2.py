# -*- coding: iso-8859-1 -*-
'''todo:
    sistemare tabella,vedi notedoc.txt,
    togliere la demenziale area@fov,mettere quella in asse -FATTO
    importando il modulo dice:
   C:\work\pythonAnalyzer\geometria\geoShell.py:1: DeprecationWarning:
   Non-ASCII character '\xb0' in file I:\vince\prog\offaxis\OA\simX2006\generaPlot.py
   on line 211, but no encoding declared;
   see http://www.python.org/peps/pep-0263.html for details 
   forse:
   possibilita' di aggiungere numero iniziale ai nomi file, in modo da averli ordinati
   routine per settare il titolo nei grafici, i valori di focale ecc, nella tabella,
   ed eventualmente il nome dei file di output
modifiche:
5-10-2006
    generaPlot2 consente di impostare resdir directory per gli output.
    aggiunto areaCopyFlag che, se vero, fa copiare i file con le aree efficaci in resdir,
    rinominandoli con il nome della configurazione .txt
istruzioni:
   per ogni cartella elencata valuta peso, campo di vista interpolato all'energia entarget
   area e angoli usati per il calcolo.
   genera file 'plotta.plt' per il plot con gnuplot di curve di area efficace e dati valutati.
   scrive gli stessi dati in un file 'tabella.txt'
   ATTENZIONE: in FOV_energie c'e' un programma come questo che oltre al FOV plotta altre
   quantita': GRASP, ecc.. i risultati vanno poi plottati con plotRoutines.
   partendo da elenco di cartelle di risultati (ottenute con i programmi fortran offAxis)
   nel di file di testo in fl,
per usarlo:
   impostare enTarget
   impostare resdir
   impostare in fl il file di testo contenente gli elenchi di directory
note:
   i valori di area efficace e peso indicato nei grafici e nelle tabelle sono rispettivamente,
   ridotto dei fattori spiderFactor (per lo spider) e aumentato del fattore
   Wfactor (per la struttura).
   '''

import os
import string
import extractInfo
from userKov.pyGeneralRoutines import  generalRoutines
import shutil
from userKov.pyOffAxis import extractInfo

if __name__=="__main__":
    #--------------------------------------------------
    # parametri
    os.chdir(r"C:\work\copiaOA\batch\batch4\batch5\simX2007")
    fileLista='lista.txt'
    entarget=30. 
    areaCopyFlag=1
    Wfactor=1.69*1.3 #fattore di moltiplicazione per i pesi in figure E tabella
    taroccoFactor=1.00 #fattore di moltiplicazione per il FOV per taroccare i dati
    spiderFactor=0.9
    resdir="FOV_energy_old"
    requirements=[] ###{0.8:100,2:1000,8:600,20:450,40:450,70:100}
    tabella=[]
    count=0   
    #---------------------------------------------------
    # inizializzazione dei files
    if not os.path.exists(resdir): os.mkdir(resdir)
    shutil.copyfile(fileLista,os.path.join(resdir,fileLista))
    plt=open(os.path.join(resdir,"plottaAll.plt"),"w")

    #se si ci sono requirements crea la stringa per aggiungerli al plot    
    reqString="\n"
    if requirements!=[]:
        rf= open(os.path.join(resdir,"requirements.txt"),"w")
        for k,v in requirements.items(): rf.write("%s\t%s\n"%(k,v))
        reqString=",\\\n 'requirements.txt' title 'Requirements' w p"
        rf.close()

    def areaFile (resdir,config,areaCopyFlag):
        ''' restituisce il nome del file da plottare, lo copia con il nome della
        configurazione se areaCopyFlag e' vero '''
        origFileName = config  + "\\aree.txt"
        destFileName = os.path.join(resdir,config + ".txt")
        if not areaCopyFlag:
            return os.path.join(os.path.pardir,origFileName)
        else:
            shutil.copyfile(origFileName,destFileName)            
            return config + ".txt"

    listaDir=open(fileLista,"r")
    a=listaDir.read().split()
    for config in a:
        print count+1, config
        Weight,ang,FOV,onAx=extractInfo.offAxisInfo(config ,entarget)[0:4]
        tabella.append([config ,Weight*Wfactor,FOV*2*taroccoFactor,onAx*spiderFactor,ang])
        plt.write("set title ' %s, ang. shell sep.=%s, FL= %s m, Dmax=%s cm' \n"  %tuple(extractInfo.extractTitle(config)))
        plt.write("set term win\n")
        plt.write("set out\n")
        plt.write("set key box\n")
        plt.write("set grid\n")
        plt.write("unset label\n")
        plt.write("set xlabel 'Energy (keV)'\n")
        plt.write("set ylabel 'Effective Area (cm^2) x 90%'\n")
        plt.write('set label "WEIGHT (+struct): %4.2f Kg\\nFOV@%4.2f keV: %4.2f arcmin" at graph 0.66,0.50\n' %(Weight*Wfactor,entarget,FOV*2*taroccoFactor))
        dataFile=areaFile(resdir,config,areaCopyFlag)
        plt.write("plot '%s' u 1:($2*%s) index 0 w l title '%s arcmin',\\"%(dataFile,spiderFactor,str(ang[0])))   
        for i in range(1,len(ang)):
            plt.write("\n'%s' u 1:($2*%s) index %s w l title '%s arcmin',\\"%(dataFile,spiderFactor,str(i),str(ang[i])))
        plt.write("\n'%s' u 1:($2*0.5*%s) index 0 w l lw 2 title '50%% on Axis Area'" %(dataFile,spiderFactor))
        plt.write(reqString)
        plt.write("\n")    
        plt.write("pause -1\n")
        plt.write("set term png\n")
        plt.write("set output '"+config +".png'\n")
        plt.write("replot\n")
        plt.write("set term win\n")
        plt.write("set out\n")
        plt.write("\n")
        count=count+1
    listaDir.close()

    plt.close()

    tab=open(os.path.join(resdir,"tabella.txt"),"w")
    tab.write("Coating\tangular shell sep (°) \tFocal lenght(m)\tMax diameter(cm)\t")
    tab.write("Weight(Kg) with struct\tFOV(arcmin diam)\tArea(cm^2 on Axis)\tSampling Angles(arcmin)\n")
    for el in tabella:
        tab.write("%s\t%s\t%s\t%s\t%2.2f\t%2.2f\t%2.2f\t%s\n"%tuple(extractInfo.extractTitle(el[0])+el[1:]))
    tab.close()
    
