# -*- coding: latin-1 -*-
'''generaPlot3 unisce generaPlot2 e FOV_energy'''
'''generaPlot:
todo:
    aggiungere nella tabella FOV/area efficace a diverse energie
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
    per ogni cartella elencata valuta peso, campo di vista interpolato all'energia enTarget
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

'''FOV_energy:
calcola fov, grasp e grasp/platescale per i file (risultanti dal programma fortran offAxis)
elencati nei file di testo in lista fl. i risultati sono nei file di testo con 4 colonne:
energia, fov diameter(arcmin), grasp (FOV/Aeff**2),grasp / plateScale(Rspot(=15")/Focal).
attenzione alla definizione della griglia di punti per l'output in enTarget.
I valori sono calcolati fino ad un energia massima per la quale si ha Aeff(Emax)=cutoffEn*Aeff[E0],
i valori di taglio vengono elencati in cutoffEnergy.txt.
Il programma crea un file per plottare le singole FOV con peso e area a 30 kev, per plottare 
gli altri risultati usare le routines in plotRoutines.py'''

import os
import string
import extractInfo
from userKov.pyGeneralRoutines import  generalRoutines
import shutil
from userKov.pyOffAxis import extractInfo



def generaPlot(rootdir,resdir,fileLista="lista.txt",enTarget=30.,requirements=[],
               spiderFactor=0.9,Wfactor=1.3,areaCopyFlag=1,taroccoFactor=1,cutoffRatio=0.05):
    '''per ogni cartella elencata valuta peso, campo di vista interpolato all'energia enTarget
    area e angoli usati per il calcolo.
    genera file 'plotta.plt' per il plot con gnuplot di curve di area efficace e dati valutati.
    scrive gli stessi dati in un file 'tabella.txt'''
    #--------------------------------------------------
    # parametri
    tabella=[]
    count=0
    tabFile="tabella_1.txt"
    #---------------------------------------------------
    # inizializzazione dei files
    os.chdir(rootdir)
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
        Weight,ang,FOV,onAx=extractInfo.offAxisInfo(config ,enTarget)[0:4]
        tabella.append([config ,Weight*Wfactor,FOV*2*taroccoFactor,onAx*spiderFactor,ang])
        plt.write("set title ' %s, ang. shell sep.=%s, FL= %s m, Dmax=%s cm' \n"  %tuple(extractInfo.extractTitle(config)))
        plt.write("set term win\n")
        plt.write("set out\n")
        plt.write("set key box\n")
        plt.write("set grid\n")
        plt.write("unset label\n")
        plt.write("set xlabel 'Energy (keV)'\n")
        plt.write("set ylabel 'Effective Area (cm^2) x 90%'\n")
        plt.write('set label "WEIGHT (+struct): %4.2f Kg\\nFOV@%4.2f keV: %4.2f arcmin" at graph 0.66,0.50\n' %(Weight*Wfactor,enTarget,FOV*2*taroccoFactor))
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

    tab=open(os.path.join(resdir,tabFile),"w")
    tab.write("Coating\tangular shell sep (�) \tFocal lenght(m)\tMax diameter(cm)\t")
    tab.write("Weight(Kg) with struct\tFOV(arcmin diam)\tArea(cm^2) on Axis\tSampling Angles(arcmin)\n")
    for el in tabella:
        tab.write("%s\t%s\t%s\t%s\t%2.3f\t%2.3f\t%2.3f\t%s\n"%tuple(extractInfo.extractTitle(el[0])+el[1:]))
    tab.close()


def calcolaFOV_Grasp(rootdir,resdir,fileLista,enTarget,requirements,spiderFactor,Wfactor,areaCopyFlag,
               taroccoFactor,cutoffRatio):
    '''calcola i risultati di FOV e GRASP per le cartelle (risultati di offaxis)
    elencate nel file fl. i risultati vanno in resdir
    cutoffRatio rapporto tra aeffMax e aeff(E) al di sotto del quale si ritiene
    non affidabile il calcolo del FOV'''
    #----------------------------------------------
    #inizializzazione
    os.chdir(rootdir)
    if not os.path.exists(resdir): os.mkdir(resdir)
    plt=open(os.path.join(resdir,"FOVvsEnergy.plt"),"w")
    cut=open(os.path.join(resdir,"cutoffEnergy.txt"),"w")
    tab=open(os.path.join(resdir,"tabella.txt"),"w")
    tabella=[]
    f=open(os.path.join(resdir,"note.txt"),"w")
    f.write(__doc__)
    f.close()
    print __doc__
    #----------------------------------------------

    listaDir=open(fileLista,"r")
    a=listaDir.read().split()
    for config in a:
        print config
        ener=extractInfo.loadEn(os.path.join(config,"aree.txt"))
        #enTarget=ener[:]
        aree=generalRoutines.loadCol(os.path.join(config,"aree.txt"),1)[0] #recupera vettore delle aree in asse
        m=max(aree)
        cutoffEn=ener[[i>m*cutoffRatio for i in aree].count(True)]
        cut.write(config+"\t"+str(cutoffEn)+"\n")
        enTargetCut=ener[0:[i<=cutoffEn for i in ener].count(True)]
        FOV_vec_diam=[2*taroccoFactor*extractInfo.offAxisInfo(config,e)[2] for e in enTargetCut]
        Weight,ang,FOV30,onAx,fovAr=extractInfo.offAxisInfo(config,enTarget)
        tabella.append([config,Weight*Wfactor,FOV30*2*taroccoFactor,onAx*spiderFactor,ang])
            
        outf=open(os.path.join(resdir,config+"_grasp.txt"),"w")
        grasp=[(FOV_vec_diam[i]**2)*generalRoutines.interpola(ener,aree,enTargetCut)[i] for i in range(len(FOV_vec_diam))]
        grasp2=[float(extractInfo.extractTitle(config)[2])*g/15 for g in grasp]
        for i in range (len(enTargetCut)):
            outf.write(str(ener[i])+"\t"+str(FOV_vec_diam[i])+"\t"+str(grasp[i])+"\t"+str(grasp2[i])+"\n")
        outf.close()
        
        plt.write("set title ' %s, ang. shell sep.=%s, FL= %s m, Dmax= %s cm' \n"  %tuple(extractInfo.extractTitle(config)))
        plt.write("set term win\n")
        plt.write("set out\n")
        plt.write("set key box\n")
        plt.write("set grid\n")
        plt.write("set xrange [0:80]\n")
        plt.write("set yrange [0:18]\n")
        plt.write("unset label\n")
        plt.write("set xlabel 'Energy (keV)'\n")
        plt.write("set ylabel 'Interpolated FOV diameter (arcmin)'\n")
        label='set label "WEIGHT (+30%%): %4.2f Kg\\n'%(Weight*Wfactor)
        label=label+ 'FOV@%4.2f keV: %4.2f arcmin\\n' %(enTarget,FOV30*2*taroccoFactor)
        label=label+ 'On-axis Aeff@%4.2f keV: %4.2f cm^2" at graph 0.02,0.15\n' %(enTarget,onAx*spiderFactor)
        plt.write(label)        
        plt.write("plot '" + config + "_grasp.txt' u 1:2 w l title '"+config+"'\n")
        
        plt.write("pause -1\n")
        plt.write("set term png\n")
        plt.write("set output 'FOVvsEn"+config+".png'\n")
        plt.write("replot\n")
        plt.write("set term win\n")
        plt.write("set out\n")
        plt.write("\n")
    listaDir.close()

    plt.close()
    cut.write("Aeff cutoff Ratio: "+str(cutoffRatio))
    cut.close()

    tab.write("Coating\tangular shell sep (�) \tFocal lenght(m)\tMax diameter(cm)\t")
    tab.write("Weight(Kg)with structure\tFOV diam(arcmin)\tArea on Axis(cm^2) x %s\tSampling Angles(arcmin)\n"%spiderFactor)
    for el in tabella:
        tab.write("%s\t%s\t%s\t%s\t%2.3f\t%2.3f\t%2.3f\t%s\n"%tuple(extractInfo.extractTitle(el[0])+el[1:]))
    tab.close()


def calc(rootdir,resdir,enTarget,fileLista,spiderFactor,Wfactor,areaCopyFlag,cutoffRatio=0.05,
         requirements=[],taroccoFactor=1.00):

    generaPlot(rootdir,resdir,fileLista,enTarget,requirements,spiderFactor,Wfactor,areaCopyFlag,
               taroccoFactor,cutoffRatio)
    calcolaFOV_Grasp(rootdir,resdir,fileLista,enTarget,requirements,spiderFactor,Wfactor,areaCopyFlag,
              taroccoFactor,cutoffRatio)
    #plottaGrasp(resdir,(0,2,3))
    #plottaFOV(resdir,(0,2,3))
    #plottaGrasp(resdir,(0,1))

if __name__=="__main__":
    
    ##simbolx
    calc(rootdir=r"H:\vince\work\gliOA\copiaOA\simX2006",
         resdir="FOV_energy_30m",
         enTarget=30.,
         fileLista='lista.txt',
         spiderFactor=0.9,
         Wfactor=1.3, #fattore di moltiplicazione per i pesi in figure e tabella (1.69*1.3 ulitmi res)
         areaCopyFlag=1,
         taroccoFactor=1.00, #fattore di moltiplicazione per il FOV per taroccare i dati
         requirements=[], ###{0.8:100,2:1000,8:600,20:450,40:450,70:100}
         cutoffRatio=0.05)
    '''
    # EDGE
    calc(rootdir=r"C:\work\copiaOA\edge_eq",
    resdir="FOV_energy",
    enTarget=2.,
    fileLista='lista.txt',
    spiderFactor=0.9,
    Wfactor=1.3, #fattore di moltiplicazione per i pesi in figure e tabella
    areaCopyFlag=1,
    taroccoFactor=1.00, #fattore di moltiplicazione per il FOV per taroccare i dati
    requirements=[], #{0.8:100,2:1000,8:600,20:450,40:450,70:100}
    cutoffRatio=0.05)'''
    '''
    # SVOM    
    calc(rootdir=r"C:\work\gliOA\copiaOA\SVOM_01",
    resdir="FOV60_svom01",
    enTarget=1.,
    fileLista='lista.txt',
    spiderFactor=1,
    Wfactor=1, #fattore di moltiplicazione per i pesi in figure e tabella
    areaCopyFlag=1,
    taroccoFactor=1.00, #fattore di moltiplicazione per il FOV per taroccare i dati
    requirements=[], #{0.8:100,2:1000,8:600,20:450,40:450,70:100}
    cutoffRatio=0.05)
    '''
    '''
    # SVOM    
    calc(rootdir=r"C:\work\gliOA\copiaOA\SVOM_01",
    resdir="FOV60_svom01",
    enTarget=1.,
    fileLista='lista.txt',
    spiderFactor=1,
    Wfactor=1, #fattore di moltiplicazione per i pesi in figure e tabella
    areaCopyFlag=1,
    taroccoFactor=1.00, #fattore di moltiplicazione per il FOV per taroccare i dati
    requirements=[], #{0.8:100,2:1000,8:600,20:450,40:450,70:100}
    cutoffRatio=0.05)
    '''
    '''
    # eROSITA   
    calc(rootdir=r"D:\2\erosita",
    resdir="FOV60_svom01",
    enTarget=1.,
    fileLista='lista.txt',
    spiderFactor=0.9,
    Wfactor=1.3, #fattore di moltiplicazione per i pesi in figure e tabella
    areaCopyFlag=1,
    taroccoFactor=1.00, #fattore di moltiplicazione per il FOV per taroccare i dati
    requirements=[], #{0.8:100,2:1000,8:600,20:450,40:450,70:100}
    cutoffRatio=0.05)      '''