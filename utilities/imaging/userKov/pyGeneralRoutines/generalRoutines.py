
import math
from numpy import interp
import copy
#from pylab import *
import os
import sys


def readFile(file,skip=0,comment=[]):
    '''return a list of strings with the lines starting from
    <skip> and with comments removed'''
    def filterCom(st,comList):
        for com in comList:
            if (st.find(com) !=-1): st=st[:st.find(com)]
        return st
     
    lines=open(file,'r').readlines()[skip:]
    lines=[filterCom(l,comment) for l in lines]
    lines=[l.strip() for l in lines]
    return lines

def readArray(file,skip=0, comment=[]):
    '''read a file in columns and load values in a matrix, skipping the
    first <skip> lines and the lines beginning with a character in 
    <comment> character'''
    lines=readFile(file)
    mat=[map(float,l.split()) for l in lines]
    return array(mat)

def locate(ar,val,sameLen=False):
    '''ar e' un vettore di numeri ordinati.
    restituisce l'indice del primo elemento maggiore di val.
    se <sameLen> e' settato restituisce l'indice dell'ultimo elemento anche se val e' maggiore,
    altrimenti restituirebbe len(ar).
    '''
    minInd=len([n for n in ar if n < val])
    assert (minInd >= 0) and (minInd<=len(ar))
    if sameLen: minInd=min(minInd,len(ar)-1)
    return minInd    

def loadCol(file,n,comment="",skip=0,nlines=None):
    ''' carica da un file l'n-sima colonna (base1), restituisce lista
    di liste di valori letti (una per ogni blocco di dati). nlines se settato
    dice quante linee leggere, se negativo quante linee dalla fine ignorare.'''
    #print 'attenzione controlla base 0 !'
    
    f=open(file,"r")
    lin=f.readlines()
    if nlines<0:
        endlin=nlines
    elif nlines == None:
        endlin=len(lin)
    else:
        endlin=skip+nlines+1
    lin=lin[skip:endlin]
    i=0
    area = []
    flag = False    #True se la precedente contiene dati, se dato dopo bianco crea nuovo blocco
    #print area
    lin=[ll for ll in lin if (ll.strip()+' ')[0]!=comment] #l'aggiunta dello spazio e' uno stupido trucco per evitare
                                                            #errore con stringa vuota
    for l in lin:
        if l.strip()=="":
            flag=False
        else:
            if not flag:
                area.append([])
                flag=True
            print "*",area, l
            area[-1].append(float(l.split()[n]))           
    f.close()
    if len(area)==1:area=area[0]
    
    return area

def arrPr(array, nel=3):
    print "len: %s\nelements:%s,...,%s"%(len(array),array[:nel],array[-nel:])

def rPars(points):
    '''restituisce m e q di una retta passante per i punti SHpoints=((x0,x1),(y0,y1))'''
    ((x0,x1),(y0,y1))=points
    m=(float(y1)-float(y0))/(float(x1)-float(x0))
    q=float(y1)-m*float(x1)
    return m,q
    
def plottaFile(file, cols=[0,1],title=None):
    '''plotta i dati contenuti nelle colonne della lista cols [x,y1,y2,...]
    nel file cols (per ora max 2 colonne)'''
    try:
        f=open(file).read().split("\n")
        f=[map(float,ff.split()) for ff in f if ff]
        f=array(f)
        plot (f[:,cols[0]],f[:,cols[1]])
        if title: legend(title)
        show ()
        return 0
    except:
        return -1
    
def searchByVal(d,val):
    '''restituisce la prima chiave trovata con valore val, restituisce None
    se non la trova'''
    for k,v in d.items():
        if v==val:return k
    return None

def sortByVal(d):
    '''restituisce due liste chiavi e valori delle chiavi e valori di d
    ordinate in base ai valori'''
    dd=d.copy()
    v=d.values()
    v.sort()
    k=[]
    for vv in v:
        a=searchByVal(d,vv)
        k.append(a)
        del d[a]     
    return k,v
    
    
def ciclazza (possibilita, indici=None, base=None):
    '''data una lista possibilita di liste di elementi tra cui scegliere, restituisce
    tutte le possibili scelte ciclando sugli indici nella tupla indici.
    per gli elementi non ciclati usa i valori in base. Senza base usa il primo, mentre
    passando possibilita come base li lascia immutati'''
    
    combs=[]
    if (indici == None): indici=range(len(possibilita))
    indici=tuple(indici)
    if base==None:
        b=[[0]*(len(possibilita))]
    else:
        b=[base]
    for i in indici:
        for bb in b:                    #print "ricavo nuovi valori da bb ",bb
            for x in possibilita[i]:
                bbb=bb[:]
                bbb[i]=x                #print "aggiungo ",bbb," a combs: ",combs
                combs.append(bbb)
            b=combs[:]
        combs=[]
    return b

'''eucDist,distanceTable e QTclust copiate dalle versioni in pythonAnalyzer '''
def eucDist(p1,p2):
    '''demenziale distanza euclidea'''
    if isinstance(p1,float):
        p1=[p1]
    if isinstance(p2,float):
        p2=[p2]
    if len(p1)<>len(p2): raise "different len p1p2"
    d=0
    for i in range(len(p1)):
        d=d+(p1[i]-p2[i])**2
    d=math.sqrt(d)
    return d

defaultDist=eucDist

def distanceTable(pList,dist=defaultDist):
    '''da una lista di punti pList, crea la tabella delle distanze con
    la distanza dist'''
    ndim=len(pList)
    t=[dist(u,v) for u in pList for v in pList]
    t=resize(t,(ndim,ndim))
    return t

def QTclust(puntiDic,Dmax,distance=defaultDist):
    '''puntiDic e' un dizionario con indice come chiave e valori su cui calcolare
    la distanza distance. Dmax e' il raggio massimo del cluster.
    restituisce una lista di liste ognuna delle quali e' un gruppo di indici che
    formano un cluster. usa algoritmo "quality threshold clustering (QTclust)"
    in caso di cluster candidati con stesso numero di membri sceglie quello di raggio
    minore'''

    def clustRad(cc,n,disTab):
        '''dato un candidato cluster cc e l'indice n di un punto nella tavola delle
        distanze disTab, restituisce la massima distanza del punto n da ogni punto del
        cluster.'''
        SelCol=choose(cc,disTab[n])
        m=reduce(maximum,SelCol)
        return m

    punti=puntiDic.copy()    
    np=len(punti)
    p=[v for v in punti.values()] #p e' una lista con i valori
    k=[v for v in punti.keys()] #k con le chiavi
    disTab=distanceTable(p,distance)
    clusterList=[]
    max_ccm=0
    minrad=0
    for j in range(np):
        cc=[] #candidate cluster
        ccrad=0
        col=disTab[j]
        interestArg= nonzero(less(col,Dmax))[0]    #array di indici dei valori che possono essere presi in considerazione
            #chissa' perche' nonzero (ma anche where) restituisce array bidimensionale...
        if len(interestArg)<max_ccm:          #e' il massimo numero teorico di membri per questo cluster
            continue     #se meno di quelli del miglior cluster candidato..
        #crea un nuovo array solo con i valori da considerare, interestArg tiene conto della corrispondenza
        #con l'array originario
        interestVal= compress(less(col,Dmax),col) #array di valori che possono essere presi in considerazione
        sv=argsort(interestVal) #indici ordinati dal piu' piccolo al + grande

        #crea cluster candidato
        cc.append(j)
        for cp in sv[1:]:
            cpArg=interestArg[cp]   #cpArg indice corrispondente nel vettore punti
            cr=clustRad(cc,cpArg,disTab)
            if cr<=Dmax:
                cc.append(interestArg[cp])
                if cr>ccrad:ccrad=cr
        ccm=len(cc)

        #decide se e' il candidato migliore finora trovato
        if ccm>max_ccm or (ccm==max_ccm and ccrad<minrad):
            max_ccm=ccm
            ccc=cc
            minrad=ccrad
            #potrebbe celare ambiguita', cluster con stesso numero di membri,ma membri diversi e
            #stesso raggio, non garantisce unicita' del risultato
    clusterList.append([k[i] for i in ccc])
    for i in ccc: del punti[k[i]]
    if len(punti)<>0: clusterList.extend(QTclust(punti,Dmax,distance))
    return clusterList
            
if __name__=='__main__':
    def test_loadCol():
        cd=os.path.dirname(sys.argv[0])
        print "-----------------------\ntest loadCol"
        testdir=os.path.join(cd,'test')
        testlist=[f for f in os.listdir(testdir) if os.path.splitext(f)[-1] == '.txt']
        for j,file in enumerate(testlist):
            f=loadCol(os.path.join(testdir,file),-1)
            print "test reading file %s: %s"%(j,file)
            print f
        print "-----------------------"
    def test_readFile():
        print readFile(r'test\readfile01.txt',0,['#'])
    def test_readArray(pwd):
        return readArray(os.path.join(pwd,'test','readArray01.txt'))

    pw=os.path.dirname(sys.argv[0])
    a=test_readArray(pw)
                         