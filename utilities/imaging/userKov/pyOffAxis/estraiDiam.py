import os
import string

def parseTransRayTrace(folder,file="transRayTrace.dat"):
    '''analizza file transraytrace.dat generato da geoinp nella cartella folder,
    per restituire lista di liste e lista di etichette. non testato con label su piu' linee,
    ne' con piu' valori numerici. non testato se il file termina senza riga bianca'''
    
    labels=[]
    values=[]
    file=os.path.join(folder,file)
    linee=open(file,"r").read().split("\n")
    status="string" #serve per capire se estendere o creare nuovo elemento in lista
    s=""
    valEl=[]
    for l in linee:
        if status=="string":
            try:
                float(l.split()[0])
                #la linea inizia con un numero
                labels.append(s.strip())
                s=""
                valEl.append(l)
                status="float"
            except:  #la linea inizia con una stringa     
                s=s+l
        elif status=="float": 
            try:
                float(l.split()[0])
                valEl.append(l)
            except:
                values.append(valEl)
                valEl=[]
                s=s+l
                status="string"
    if s != "":labels.append(s.strip())
    if valEl!=[]:values.append(valEl)
    return values, labels
                
def estraiDiamDic(dirList,basedir=os.getcwd()):
    '''passando una lista di nomi di directory contenute in basedir
    e ognuna contenente i file di un risultato di transRayTrace,
    restituisce un dizionario che ha per chiavi i percorsi dei file e
    per valori le liste dei diametri ATTENZIONE CANNA DI 1 shell'''
    tmp=os.getcwd()
    os.chdir(basedir)
    diamDic={}
    for d in dirList:
        diamList=[s.strip() for s in parseTransRayTrace(d)[0][3]]
        diamDic[d]=diamList
    os.chdir(tmp)
    return diamDic    

def estraiDiam( fileLista="estraiDiam.txt",invertiDiam=False):
    '''estrae i diametri dai file transRayTrace.dat per directory elencate
    nel file fileLista. li mette nella cartella outDiam (nella dir di fileLista)
    con il nome del file uguale a quello della directory originaria '''

    #prepara le directory
    baseDir=os.path.split(fileLista)[0]
    flName=os.path.split(fileLista)[1]  
    if baseDir == '': baseDir=os.path.curdir
    outdir=os.path.join(baseDir,"extractedDiam")
    dl=open(fileLista,"r").read().split()
    count=0
    diamDic=estraiDiamDic(dl,baseDir)
    if not (os.path.isdir(outdir)):
        os.mkdir(outdir)  
    for d in dl:
        diam=diamDic[d]
        c_out=open(os.path.join(outdir,d+".txt"),"w")
        if invertiDiam: diam=diam[::-1]
        cc="\n".join(diam)
        c_out.write(cc)
        c_out.close()
        count=count+1

    return "estratti diametri da %s cartelle\nrisultati in %s" %(count,outdir)

def mostraLimiti(fileLista="estraiDiam.txt"):
    '''restituisce stringa con output: diametro min e max per ogni configurazione
    in fileLista. fileLista deve trovarsi nella directory in cui si trovano
    le cartelle dei risultati di offAxis elencate'''
    
    baseDir=os.path.split(fileLista)[0]
    fileName=os.path.split(fileLista)[1]
    if baseDir == '': baseDir=os.path.curdir
    dl=open(fileLista,"r").read().split()
    dl=[d for d in dl if d[-2:]=="60"]
    d=estraiDiamDic(dl,baseDir)
    sl=["folder\tdmin\tdmax"]
    a=d.keys()
    a.sort()
    for k in a:
      sl.append("%s\t%s\t%s"%(k,d[k][-1],d[k][0]))
    return "\n".join(sl)
    
    
if __name__=="estraiDiam":
    print "questo e' estraiDiam"

if __name__=="__main__":
    print mostraLimiti(r"C:\work\copiaOA\batch\simx2007\lista.txt")


            
                

