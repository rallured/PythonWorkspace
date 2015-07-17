# -*- coding: latin-1 -*-
'''
4.4 non cambiato niente, volevo farlo funzionare con fovfraction vettoriale.
22-11-2008
generaplot4.3
changed the data in <filename>_grasp.txt. Now 1st col is fov diam at FOVfraction, 2nd is
FOV at FOVfraction2, 3rd is grasp in cm^2 deg^2, 4th is platescale.
The new definition of grasp is 0.75*onaxisaeff*pi*fovradius^2
17-11-2008
added a variable fraction area, to indicate the fraction of area to
be assumed for the FOV value (typical 0.5 for FOV, 0.8 for GRASP)
13-6-2008 
generaPlot4.2:
--!! The results from this last version (tested on  H:\vince\work\gliOA\copiaOA\simX2006)
seems to be different (and more correct!) from the previous ones in that folder. The FOV
vs energy doesn't change, but the values in the labels were wrong.
-added information in log file, like time,date version and parameters.
-added a reminder to update the version.
-replaced the 'interpola' routine with numpy.interp
-added 'with spider' in the fov vs energy label on plots
27-4-07 tentativo non concluso di aggiungere nella lista
una seconda colonna opzionale con la cartella da cui recuperare
aree in asse, in modo da poter fare grafici anche per intervalli
di angoli off Axis non comprendenti lo zero.

13-4-07 generaPlot4: modificato il main con passaggio parametri a routine
modificato il conteggio di enercut e taglio dei vettori, non dovrebbe dare
piu' errori se include tutti i valori ne' se termina prima.

generaPlot3 unisce generaPlot2 e FOV_energy
'''

'''
generaPlot:
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
import math
from userKov.pyGeneralRoutines import  generalRoutines
import shutil
from userKov.pyOffAxis import extractInfo
from userKov.pyOffAxis import plotRoutines
from time import strftime
from numpy import interp

version = "generaPlot4.4"

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
    a=listaDir.read().split("\n")
    for config in a:
        if len(config.split())==1:
            onAxisFolder=None
        else:
            config,onAxisFolder=config.split()
        print count+1, config
        Weight,ang,FOV,onAx=extractInfo.offAxisInfo(config ,enTarget,onAxisFolder,fractionArea=0.5)[0:4]
        tabella.append([config ,Weight*Wfactor,FOV*2*taroccoFactor,onAx*spiderFactor,ang])
        plt.write("set title ' %s, ang. shell sep. = %s�, FL = %s m, Dmax = %s cm' \n"  %tuple(extractInfo.extractTitle(config)))
        plt.write("set term win\n")
        plt.write("set out\n")
        plt.write("set key box\n")
        plt.write("set grid\n")
        plt.write("unset label\n")
        plt.write("set xlabel 'Energy (keV)'\n")
        plt.write("set ylabel 'Effective Area (cm^2) x 90%'\n")
        plt.write('set label "WEIGHT (+struct): %4.2f kg\\nFOV@%4.2f keV: %4.2f arcmin" at graph 0.66,0.50\n' %(Weight*Wfactor,enTarget,FOV*2*taroccoFactor))
        dataFile=areaFile(resdir,config,areaCopyFlag)
        plt.write("plot '%s' u 1:($2*%s) index 0 w l lw 2 title '%s arcmin',\\"%(dataFile,spiderFactor,str(ang[0])))   
        for i in range(1,len(ang)):
            plt.write("\n'%s' u 1:($2*%s) index %s w l lw 2 title '%s arcmin',\\"%(dataFile,spiderFactor,str(i),str(ang[i])))
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
    tab.write("Weight(kg) with struct\tFOV(arcmin diam)\tArea(cm^2) on Axis\tSampling Angles(arcmin)\n")
    for el in tabella:
        tab.write("%s\t%s\t%s\t%s\t%2.3f\t%2.3f\t%2.3f\t%s\n"%tuple(extractInfo.extractTitle(el[0])+el[1:]))
    tab.close()


def calcolaFOV_Grasp(rootdir,resdir,fileLista,enTarget,requirements,spiderFactor,Wfactor,areaCopyFlag,
               taroccoFactor,cutoffRatio,FOVfraction=0.5,FOVfraction2=0.8):
    '''calcola i risultati di FOV e GRASP per le cartelle (risultati di offaxis)
    elencate nel file fl. i risultati vanno in resdir.
    spiderFactor non incide sui dati numerici, ma solo sui plot e sui valori in tabella. 
    cutoffRatio rapporto tra aeffMax e aeff(E) al di sotto del quale si ritiene
    non affidabile il calcolo del FOV'''
    #----------------------------------------------
    #inizializzazione
    os.chdir(rootdir)
    plsize=15
    if not os.path.exists(resdir): os.mkdir(resdir)
    plt=open(os.path.join(resdir,"FOVvsEnergy.plt"),"w")
    cut=open(os.path.join(resdir,"cutoffEnergy.txt"),"w")
    tab=open(os.path.join(resdir,"tabella.txt"),"w")
    tabella=[]
    f=open(os.path.join(resdir,"note.txt"),"w")
    f.write(__doc__)
    f.write("time (Y/M/D): "+strftime("%Y-%m-%d %H:%M:%S"))
    f.write("""-------\nrootdir: %s\nresdir: %s\nfileLista: %s\nenTarget: %s\nrequirements: %s\n
            spiderFactor: %s\nWfactor: %s\nareaCopyFlag: %s\ntaroccoFactor: %s\ncutoffRatio: %s"""
            %(rootdir,resdir,fileLista,enTarget,requirements,spiderFactor,Wfactor,areaCopyFlag,
            taroccoFactor,cutoffRatio))
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
        cutIndex=[i>m*cutoffRatio for i in aree].count(True)
        enTargetCut=ener[:cutIndex]
        cutoffEn=enTargetCut[-1]
        cut.write(config+"\t"+str(cutoffEn)+"\n")
        FOV_vec_diam=[2*taroccoFactor*extractInfo.offAxisInfo(config,e,fractionArea=FOVfraction)[2] for e in enTargetCut]
        Weight,ang,FOV30,onAx,fovAr=extractInfo.offAxisInfo(config,enTarget,fractionArea=FOVfraction)
        tabella.append([config,Weight*Wfactor,FOV30*2*taroccoFactor,onAx*spiderFactor,ang])
            
        outf=open(os.path.join(resdir,config+"_grasp.txt"),"w")
        FOV2_vec_diam=[2*taroccoFactor*extractInfo.offAxisInfo(config,e,fractionArea=FOVfraction2)[2] for e in enTargetCut]
        
        grasp=[(0.75*math.pi*(FOV_vec_diam[i]**2)/4)*interp(enTargetCut,ener,aree)[i]/3600. for i in range(len(FOV_vec_diam))]
        focale=float(extractInfo.extractTitle(config)[2])
        plscale=[focale*g/plsize for g in grasp]
        outf.write("Energy(keV)\tFOV at %s%%\tFOV at %s%% onaxis Aeff\tgrasp(0.75*aeffonaxis*FOVarea)\tplatescale(grasp*F(=%sm)/%scm)\n"%(str(FOVfraction),
                                                                                                str(FOVfraction2),
                                                                                                str(focale),
                                                                                                str(plsize)))
        for i in range (len(enTargetCut)):
            outf.write(str(ener[i])+"\t"+str(FOV_vec_diam[i])+"\t"+str(FOV2_vec_diam[i])+"\t"+str(grasp[i])+"\t"+str(plscale[i])+"\n")
        outf.close()
        
        plt.write("set title ' %s, ang. shell sep. = %s�, FL = %s m, Dmax = %s cm' \n"  %tuple(extractInfo.extractTitle(config)))
        plt.write("set term win\n")
        plt.write("set out\n")
        plt.write("set key box\n")
        plt.write("set grid\n")
        plt.write("set xrange [0:*]\n")
        plt.write("set yrange [0:*]\n")
        plt.write("unset label\n")
        plt.write("set xlabel 'Energy (keV)'\n")
        plt.write("set ylabel 'Interpolated FOV diameter (arcmin)'\n")
        label='set label "WEIGHT (+30%%): %4.2f kg\\n'%(Weight*Wfactor)
        label=label+ 'FOV@%4.2f keV: %4.2f arcmin\\n' %(enTarget,FOV30*2*taroccoFactor)
        label=label+ 'On-axis Aeff (with spider)@%4.2f keV: %4.2f cm^2" at graph 0.02,0.15\n' %(enTarget,onAx*spiderFactor)
        plt.write(label)
        ss="plot '" + config + "_grasp.txt' u 1:2 w l lw 2 title '"+str(FOVfraction)+" on axis Aeff"
        if (FOVfraction == FOVfraction2):
            ss=ss+"'\n"
        else:
            ss=ss+"',\\\n'" + config + "_grasp.txt' u 1:3 w l lw 2 title '"+str(FOVfraction2)+" on axis Aeff'\n"
        plt.write(ss)
        
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
    tab.write("Weight(kg)with structure\tFOV diam for %s %%on axis area(arcmin)\t\
              Area on Axis(cm^2)x %s\tSampling Angles(arcmin)\n" %(FOVfraction,spiderFactor))
    for el in tabella:
        tab.write("%s\t%s\t%s\t%s\t%2.3f\t%2.3f\t%2.3f\t%s\n"%tuple(extractInfo.extractTitle(el[0])+el[1:]))
    tab.close()


def calc(rootdir,resdir,enTarget,fileLista,spiderFactor,Wfactor,areaCopyFlag,cutoffRatio=0.05,
         requirements=[],taroccoFactor=1.00,onAxisFolder=None,FOVfraction=None,FOVfraction2=None):

    generaPlot(rootdir,resdir,fileLista,enTarget,requirements,spiderFactor,Wfactor,areaCopyFlag,
               taroccoFactor,cutoffRatio)
    calcolaFOV_Grasp(rootdir,resdir,fileLista,enTarget,requirements,spiderFactor,Wfactor,areaCopyFlag,
              taroccoFactor,cutoffRatio,FOVfraction,FOVfraction2)
    plotRoutines.cmpPlot(resdir,overwrite=1,groupFile=fileLista,xrange="[*:*]",yrange="[*:*]",pltFilePrefix="Grasppng_",
        outPrefix="Grasp_",term="ps",
        using="u 1:4",ylabel="Grasp (Aeff x area FOV) [cm^2 deg^2]",title="GRASP")


if __name__=="__main__":
    '''
    ##simbolx
    calc(rootdir=r"C:\work\copiaOA\batch\batch4\simx2007_100kev",
         resdir="FOV_energy",
         enTarget=30.,
         fileLista='lista.txt',
         spiderFactor=0.9,
         Wfactor=1.69*1.3, #fattore di moltiplicazione per i pesi in figure e tabella
         areaCopyFlag=1,
         taroccoFactor=1.00, #fattore di moltiplicazione per il FOV per taroccare i dati
         requirements=[], ###{0.8:100,2:1000,8:600,20:450,40:450,70:100}
         cutoffRatio=0.05)
    '''
    '''
    # EDGE
    calc(rootdir=r"C:\work\copiaOA\edge_eq3",
    resdir="FOV_energy_extr",
    enTarget=1.,
    fileLista='lista.txt',
    spiderFactor=0.9,
    Wfactor=1.3, #fattore di moltiplicazione per i pesi in figure e tabella
    areaCopyFlag=1,
    taroccoFactor=1.00, #fattore di moltiplicazione per il FOV per taroccare i dati
    requirements=[], #{0.8:100,2:1000,8:600,20:450,40:450,70:100}
    cutoffRatio=0.05)

    '''
    '''
    # hxmt
    calc(rootdir=r"C:\work\offAxis_hxmt\hxmt_12",
    resdir="FOV_energy",
    enTarget=1.,
    fileLista='lista.txt',
    spiderFactor=0.9,
    Wfactor=1.3, #fattore di moltiplicazione per i pesi in figure e tabella
    areaCopyFlag=1,
    taroccoFactor=1.00, #fattore di moltiplicazione per il FOV per taroccare i dati
    requirements=[], #{0.8:100,2:1000,8:600,20:450,40:450,70:100}
    cutoffRatio=0.05)
    '''

    '''
    calc(rootdir=r"F:\vince\risultati\EDGE\2_ray_tracing\edge_eq3",
        resdir="FOV_energy2009",
        enTarget=3.,
        fileLista='lista.txt',
        spiderFactor=0.9,
        Wfactor=1.3,        #fattore di moltiplicazione per i pesi in figure e tabella (1.69*1.3 ulitmi res)
        areaCopyFlag=1,
        taroccoFactor=1.00, #fattore di moltiplicazione per il FOV per  riscalare i dati
        ##requirements={0.8:100,2:1000,8:600,20:450,40:450,70:100}, ###requirements={0.8:100,2:1000,8:600,20:450,40:450,70:100}
        cutoffRatio=0.05,FOVfraction=0.5,FOVfraction2=0.75)'''

    '''calc(rootdir=r"F:\next_HXT",
    resdir="FOV_30keV",
    enTarget=30.,
    fileLista='fovList2.txt',
    spiderFactor=0.9,
    Wfactor=1.3,        #fattore di moltiplicazione per i pesi in figure e tabella (1.69*1.3 ulitmi res)
    areaCopyFlag=1,
    taroccoFactor=1.00, #fattore di moltiplicazione per il FOV per  riscalare i dati
    requirements={10:542,30:169,60:72}, ###requirements={0.8:100,2:1000,8:600,20:450,40:450,70:100}
    cutoffRatio=0.05,FOVfraction=0.5,FOVfraction2=0.75)'''

    '''calc(rootdir=r"F:\next_SXT",
    resdir="FOV_SXT",
    enTarget=7.,
    fileLista='listaFOV.txt',
    spiderFactor=0.9,
    Wfactor=1.3,        #fattore di moltiplicazione per i pesi in figure e tabella (1.69*1.3 ulitmi res)
    areaCopyFlag=1,
    taroccoFactor=1.00, #fattore di moltiplicazione per il FOV per  riscalare i dati
    requirements={1:611,7:469,10:305}, ###requirements={0.8:100,2:1000,8:600,20:450,40:450,70:100}
    cutoffRatio=0.05,FOVfraction=0.5,FOVfraction2=0.75)'''

    calc(rootdir=r"E:\work\workOA\traie7\hexitSat_2009",
    resdir="FOV_Energy_test",
    enTarget=30.,
    fileLista='listaFOV.txt',
    spiderFactor=0.9,
    Wfactor=1.3,        #fattore di moltiplicazione per i pesi in figure e tabella (1.69*1.3 ulitmi res)
    areaCopyFlag=1,
    taroccoFactor=1.00, #fattore di moltiplicazione per il FOV per  riscalare i dati
    #requirements={1:611,7:469,10:305}, ###requirements={0.8:100,2:1000,8:600,20:450,40:450,70:100}
    cutoffRatio=0.05,FOVfraction=0.5,FOVfraction2=0.75)
    
