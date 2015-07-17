import extractInfo
import os

def createDic(ener,aree,entarget,ang):
    '''data una lista di liste(una per ogni angolo) delle aree <aree>,
    con le aree in funzione dell'energia <ener>, restituisce un dizionario
    con chiavi gli angoli della lista <ang> degli angoli del ray-tracing, e per
    valori le aree interpolate    '''
    
    def findEnTargetIndex(ener,target):
        '''dato un vettore float delle energie ener, restituisce indice di
        quella immediatamente inferiore a target'''
        l= [i for i in ener if i < target]
        ll=len(l)
        return ll-1

    ind= findEnTargetIndex(ener,entarget)
    d={}
    for i in range(len(aree)):
        d[ang[i]]=aree[i][ind]+((aree[i][ind+1]-aree[i][ind])/(ener[ind+1]-ener[ind]))*(entarget-ener[ind])
    return d

def interpolateFOV(y,fractionArea=0.5):
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

def findFactor1(ener,aree,folder,targetFOV,angInf):
    '''determina il fattore di incremento per l'area efficace a 6 arcmin
    per avere il campo di vista targetFOV, vale per incrementi dell'angolo di 2 arcmin'''
    dic=createDic(ener,aree,30.,extractInfo.retrieveAngleSteps(folder))
    A=dict([(u,v/dic[0]) for u,v in dic.items()])
    Ainf=A[angInf]
    Asup=A[angInf+2.]
    fac=(Ainf+(1-2*Ainf)/(targetFOV/2-angInf))/Asup
    #print A
    return fac

def findFactor2(ener,aree,folder,targetFOV,angInf):
    dic=createDic(ener,aree,30.,extractInfo.retrieveAngleSteps(folder))
    A=dict([(u,v/dic[0]) for u,v in dic.items()])
    Ainf=A[angInf]
    Asup=A[angInf+2.]
    dA=Asup-Ainf
    #print A
    fac=1/(2*Ainf+dA*(targetFOV/2-angInf))    
    return fac
    
if __name__=="__main__":

    folder=r"C:\work\copiaOA\batch\batch4\simx2007_80kev"
    fileList="targetFOV.txt"

    for l in open(os.path.join(folder,fileList),"r").read().split("\n"):
        l=l.split()
        f=os.path.join(folder,l[0])
        targetFOV=float(l[1])
        aree=extractInfo.loadCol(os.path.join(f,"aree.txt"),1)
        ener=extractInfo.loadEn(os.path.join(f,"aree.txt"))
        angInf=2*float(int(targetFOV/4))
        fac1=findFactor1(ener,aree,f,targetFOV,angInf)
        fac2=findFactor2(ener,aree,f,targetFOV,angInf)
        print l[0],l[1],fac1,fac2,angInf

        kk=[0,2.0,4.0,6.0,8.0,10.0].index(angInf)
        out=open(os.path.join(folder,l[0]+"_fac1.txt"),"w")
        a=aree[:]
        for k in range(kk+1,len(a)):
            a[k]=[area*fac1 for area in a[k]]
        for aa in a:
            for t in zip(ener,aa):out.write(str(t[0])+"\t"+str(t[1])+"\n")
            out.write("\n"*2)
        out.close()

        out=open(os.path.join(folder,l[0]+"_fac2.txt"),"w")
        a=aree[:]
        for k in range(kk,len(a)):
            a[k]=[area*fac2 for area in a[k]]
        for aa in a:
            for t in zip(ener,aa):out.write(str(t[0])+"\t"+str(t[1])+"\n")
            out.write("\n"*2)        
        out.close()

        

def ricalcolaFOV(fileAree):
    '''per verificare aree risultanti'''
    #crea Dizionario normalizzato
    aree=extractInfo.loadCol(fileAree,1)
    ener=extractInfo.loadEn(fileAree)
    dic= createDic(ener,aree,30.,[0,2.0,4.0,6.0,8.0,10.0])
    A=dict([(u,v/dic[0]) for u,v in dic.items()])
    #interpolateFOV
    print fileAree, interpolateFOV(A)[0]*2, interpolateFOV(A)[1]
    print A
    print dic
        

        
    
    



    