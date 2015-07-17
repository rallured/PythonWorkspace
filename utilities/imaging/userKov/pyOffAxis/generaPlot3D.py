'''
-copiato da generaPlot per generare file con tre colonne richiestimi da malaguti

todo:
   adattare con if __main__
   forse:
   possibilita' di aggiungere numero iniziale ai nomi file, in modo da averli ordinati
   routine per settare il titolo nei grafici, i valori di focale ecc, nella tabella,
   ed eventualmente il nome dei file di output
istruzioni:
   partendo da elenchi di cartelle di risultati (ottenute con i programmi fortran offAxis)
   elencati in una lista fl di file di testo,
   per ogni cartella elencata valuta peso, campo di vista interpolato all'energia entarget
   area e angoli usati per il calcolo.
   genera file 'plotta.plt' per il plot con gnuplot di curve di area efficace e dati valutati.
   scrive gli stessi dati in un file 'tabella.txt'
per usarlo:
   impostare enTarget
   impostare in fl l'elenco dei file di testo contenenti gli elenchi di directory
   attenzione alla routine extractTitle che estrae i valori da usare nei titoli dei grafici
   e nelle tabelle dal nome delle directory (nettamente migliorabile, es. potrebbe prenderli
   dal file dei parametri)'''

import os
import string

def offAxisInfo(folder,entarget):
    '''passando un folder e l'energia entarget restituisce
    peso,lista degli angoli, fov interpolato e area corrispondente'''
    def Wrecover(folder):
        '''recupera il peso dal file pesi.txt'''
        f=open(folder+"\pesi.txt","r")
        testo=f.readlines()
        peso=testo[0].split()
        f.close()
        return float(peso[2])

    def retrieveAngleSteps(folder):
        '''recupera dal file delle impostazioni in folder i valori degli angoli offaxis utilizzati
        li restituisce in lista'''
        f=open(os.path.join(folder,"imp_OffAxis.txt"),"r")
        for l in f.readlines():
            if string.find(string.lower(l),"ang0arcmin")<>-1:
                angmin=float(l.split("=")[1])
            elif string.find(string.lower(l),"ang1arcmin")<>-1:
                angmax=float(l.split("=")[1])
            elif string.find(string.lower(l),"pasa")<>-1:
                pasa=float(l.split("=")[1])    
        f.close()
        a=angmin
        angList=[]
        while a<=angmax:
            angList.append(a)
            a=a+pasa
        return angList
        
    def evaluateFOV(folder,entarget,ang):
        '''valuta il campo di vista all'energia energy per il file delle aree
        offaxis contenuto in folder'''

        def loadEn(file):
            ''' carica da un file delle aree un vettore con le energie '''
            f=open(file,"r")
            lin=f.readlines()
            ener = []
            for l in lin:
                if l==" \n": break
                else: ener.append(float(l.split()[0]))            
            f.close()
            return ener

        def loadCol(file,n):
            ''' carica da un file l'n-sima colonna in lista di liste'''
            f=open(file,"r")
            lin=f.readlines()
            i=0
            area = [[]]
            flag = False
            #print area
            for l in lin:
                #print i,area[i]
                if l==" \n":
                    if flag:
                        i=i+1
                        flag=False
                    else:
                        flag=True
                        area.append([])
                else: area[i].append(float(l.split()[n-1]))            
            f.close()
            area.pop()
            return area
        
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


        def interpolateFOV(y):
            '''dal dizionario {x:y} con x angoli offaxis e y aree normalizzate restituisce l'angolo interpolato
            con area=fractionArea '''
            #print oaDic
            fractionArea=0.5
            k=y.keys()
            k.sort()
            for x in k:
                #print y[x]
                if y[x]>fractionArea:
                    minAng=(x,y[x])
                elif y[x]<fractionArea:
                    maxAng=(x,y[x])
                    break
                elif y[u]==fractionArea:
                    minAng=maxAng=(x,y[x])
            #print minAng,maxAng,"<--"
            fov=(((maxAng[0]-minAng[0])/(maxAng[1]-minAng[1]))*(fractionArea-maxAng[1]))+maxAng[0]
            return fov,fractionArea
        

        areeFile=os.path.join(folder,"aree.txt")
        ener=loadEn(areeFile)    
        #entarget contiene il valore di energia immediatamente inferiore a energy
        #ora serve un pezzo che prende in dizionario angolo:area interpolata
        aree=loadCol(areeFile,2)
        offAxisAr={}
        offAxisAr=createDic(ener,aree,entarget,ang)
        onAx=offAxisAr[0]
        norm=[(u,v/onAx) for u,v in offAxisAr.items()]
        norm=dict(norm)
        #print norm
        FOV,intAr = interpolateFOV(norm)
        return FOV,intAr*onAx
    
    print folder
    ang=retrieveAngleSteps(folder)
    Weight=Wrecover(folder)
    fov,fovAr=evaluateFOV(folder,entarget,ang)
    return Weight,ang,fov,fovAr


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



plt=open("plotta3d.plt","w")
file3d=open("3d.txt","w")
fl=['test3d.txt']  

for filelista in fl:
    listaDir=open(filelista,"r")
    a=listaDir.readlines()
    for config in a:
        nomeFile=os.path.join(config[:-1],"aree.txt")
        foc=extractTitle(config[:-1])[2]
        aeffFile=open(nomeFile,"r")
        vettore=aeffFile.readlines()
        vettore=vettore[0:vettore.index(' \n')]
        vettore=[v.split()[0]+"\t"+str(foc)+"\t"+v.split()[1]+"\n" for v in vettore]
        aeffFile.close()
        for v in vettore:
            file3d.write(v)
        file3d.write("\n\n")
    listaDir.close()
plt.close()

file3d.close()
