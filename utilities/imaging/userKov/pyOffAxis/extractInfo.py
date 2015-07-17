from userKov.pyGeneralRoutines import  generalRoutines
import os
import string
import namelist
'''attenzione alla routine extractTitle che estrae i valori da usare nei titoli dei grafici
e nelle tabelle dal nome delle directory (nettamente migliorabile, es. potrebbe prenderli
dal file dei parametri)'''
   
def retrieveWeights(folder):
    '''recupera il peso dal file pesi.txt in folder'''
    f=open(folder+"\pesi.txt","r")
    testo=f.readlines()
    peso=testo[0].split()
    f.close()
    return float(peso[2])

def Wrecover(folder):
    '''recupera il peso dal file pesi.txt in folder'''
    f=open(folder+"\pesi.txt","r")
    testo=f.readlines()
    peso=testo[0].split()
    f.close()
    return float(peso[2])

def extractTitle(st):
    '''dalla stringa con il nome della directory estrae le informazioni
    sulla configurazione e restituisce una lista con
    coating, filling factor, focale, diametro massimo'''
    st=str(st)
    coating={"ir":"Ir monolayer","ml":"W/Si multilayer","m2":"Pt/C multilayer",
             "m3":"Pt/C+C / Ni/C"}
    ff={"ff000":"0.0","ff00":"0.0","ff007":"0.07","ff015":"0.15"}
    diam={"max60":"60","max65":"65","max70":"70"}
    focale={"20m":"20","22m":"22.5","25m":"25","27m":"27.5","30m":"30"
            ,"180":"18.0","185":"18.5","190":"19.0","195":"19.5","200":"20.0"}
    l=[string.lower(w) for w in st.split("_")]
    try:
        k=[coating[l[0]],ff[l[1]],focale[l[2]],diam[l[3]]]
    except:
        a=namelist.read(os.path.join(st,"imp_OffAxis.txt"))
        return ""
        k=[os.path.split(st)[-1],a['CREATEDIAMS']['FieldOfViewDeg'],
            a['GEOPARS']['F_LENGTHdaImp_m'],           
            a['DIAMMAX']['massimo']]
        
    return k

def loadEn(file):
    ''' carica da un file delle aree un vettore con le energie '''
    f=open(file,"r")
    lin=f.readlines()
    ener = []
    for l in lin:
        if l.strip()=="":
            break            
        else:
            ener.append(float(l.split()[0]))            
    f.close()
    return ener

def shellAngles(folder):
    '''recupera dal file shellstruct la sequenza di angoli per le shell'''
    lista= open(os.path.join(folder,"shellStruct.txt"),"r").read().split("\n")[1:]
    ang=[riga.split()[5] for riga in lista if riga != '']
    return ang

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

def evaluateFOV(folder,entarget,ang,onAxisArea=None,fractionArea=0.5):
    '''valuta il campo di vista (per cui la frazione di area e' <fractionArea>)
    all'energia <entarget> per il file delle area offaxis in <folder>
    contenente le aree per gli angoli fuori asse della lista <ang>'''
    
    def createDic(ener,aree,entarget,ang):
        '''data una lista di liste(una per ogni angolo) delle aree <aree>,
        con le aree in funzione dell'energia <ener>, restituisce un dizionario
        con chiavi gli angoli della lista <ang> degli angoli del ray-tracing, e per
        valori le aree interpolate
        
        '''
        def findEnTargetIndex(ener,target):
            '''dato un vettore float delle energie ener, restituisce indice di
            quella immediatamente inferiore a target'''
            l= [i for i in ener if i < target]
            ll=len(l)
            return ll-1

        ind= findEnTargetIndex(ener,entarget)
        d={}
        for i in range(len(aree)):
            #print i
            d[ang[i]]=aree[i][ind]+((aree[i][ind+1]-aree[i][ind])/(ener[ind+1]-ener[ind]))*(entarget-ener[ind])
        return d
    
    def interpolateFOV(y,fractionArea):
        '''dal dizionario <y> {k:v} con k angoli offaxis e v aree normalizzate
        restituisce l'angolo interpolato con area=fractionArea '''
        #print oaDic
        k=y.keys()
        k.sort()
        maxAng=-1
        for x in k:
            #print y[x]
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
        #print minAng,maxAng,"<--"
        if maxAng[0]==minAng[0]: fov = maxAng[0]
        else:fov=(((maxAng[0]-minAng[0])/(maxAng[1]-minAng[1]))*(fractionArea-maxAng[1]))+maxAng[0]
        return fov,fractionArea
    

    areeFile=os.path.join(folder,"aree.txt")
    ener=loadEn(areeFile)    
    #entarget contiene il valore di energia immediatamente inferiore a energy
    #ora serve un pezzo che prende in dizionario angolo:area interpolata
    aree=generalRoutines.loadCol(areeFile,1)
    offAxisAr={}
    if onAxisArea:
        print "extractInfo.py: adding onAxisArea to area Dic!!"
        aree.insert(0,onAxisArea)
        ang.insert(0,0)
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
    FOV,intAr = interpolateFOV(norm,fractionArea)
    return FOV,onAx,intAr*onAx

def offAxisInfo(folder,entarget,onAxisFolder=None,fractionArea=0.5):
    '''passando un folder e l'energia entarget restituisce
    peso,lista degli angoli, fov interpolato (per cui la frazione di area e' <fractionArea>)
    e area corrispondente. onAxisFolder puo' indicare un'altra cartella contenente l'area in asse '''  
    #print folder
    if onAxisFolder:
        af=os.path.join(onAxisFolder,"aree.txt")
        print "offAxisInfo: Loading on Axis area from %s"%af
        onAxisArea=generalRoutines.loadCol(af,1)[0]
    else:
        onAxisArea=None
    ang=retrieveAngleSteps(folder)
    Weight=retrieveWeights(folder)
    fov,onAx,fovAr=evaluateFOV(folder,entarget,ang,onAxisArea,fractionArea)
    return Weight,ang,fov,onAx,fovAr


