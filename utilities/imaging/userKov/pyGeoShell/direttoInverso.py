import math

'''ricava il dianetri allintersezione dal diametro d'ingresso
con metodo ricorsivo. lunghezze con stessa unita' di misura.
verifica che sia giusto con test diretto, ricavando diametro
all' ingersso da diam allintersezione'''

def diretto (F,dmin,SH,wolter=1):
    '''ricavando diametro all' ingersso da diam all'intersezione su superficie mirror'''
    F=float(F)
    dmin=float(dmin)
    th=math.atan2(dmin/2,F)/4
    if wolter == 0: dmax=2*SH*math.tan(th)+dmin
    else : dmax=2*math.sqrt((dmin/2)**2+2*math.tan(th)*(dmin/2)*SH)
    return dmax

def inverso(F,dmax,H,tol=10e-8):
    '''ricava il diametro all'intersezione dal diametro d'ingresso
    con metodo ricorsivo. stessa unita' di misura.'''
    rIn=float(dmax)/2
    def ThSh(R):return math.atan2(float(R),float(F))/4
    def Rint(rIn,th0):return rIn-H*math.tan(th) #trova raggio all'intersezione da angolo e rIn
    r=rIn
    rOld=0
    while (abs(r-rOld)>=tol):
        rOld=r
        th=ThSh(r)        
        r=Rint(rIn,th)
    return r*2   

def nextInner(D,theta,Hpars):
    '''dati diametro all'intersezione di una shell, angolo di filling factor (ingradi)
    m e q dell'altezza shell (ottica a farfalla), t.c. SH= m*D+q,
    restituisce il diametro esterno alla pupilla d'ingresso per la shell successiva'''
    
    r0=float(D)/2
    theta=(math.pi/180)*(theta/60.)
    mH,qH=map(float,Hpars)
    rNew=(r0-qH*math.tan(theta))/(mH*math.tan(theta)+1)
    return rNew*2

def Nouter(N,D,F,Hpars,tpars,ffdeg=0,wolter=1):
    '''give the external diameter (wall included) for the Nth shell,
    outer to the given one.'''
    r0=float(D)/2
    mH,qH=map(float,Hpars)
    mt,qt=map(float,tpars)
    rNext=r0
    SH=0
    for i in range(1,N+1):
        SH=mH*rNext+qH
        rNext=diretto(F,rNext*2,SH,wolter)/2.
        t=mt*rNext*2+qt
        rNext=rNext+ t+ SH*math.tan(ffdeg*math.pi/180.) 
    return rNext*2-2*SH*math.tan(ffdeg*math.pi/180.) 

if __name__=="x__main__":
    '''ricava il dianetri all'intersezione dal diametro d'ingresso con metodo ricorsivo.
    lunghezze con stessa unita' di misura. verifica che sia giusto con test diretto,
    ricavando diametro all' ingersso da diam allintersezione'''
    
    testShell=[(3000,700,(0,70),30),(2000,700,(0,70),30),(3000, 364.964713 ,(0,100),30)]
    for (F,dmin,H,ff) in testShell:
        dmax=diretto (F,dmin,H[1])
        d=inverso(F,dmax,H[1])
        print "\n---- F=%s, SH=%s -----"%(F,H[1])
        print "dmin:%s -> dmax:%s ->\n%s"%(dmin,dmax,d)
        print "--\nnext shell \(ff=%s, H = %s*D+%s):"%(ff,H[0],H[1])
        dNext=nextInner(d,ff,H)
        print "%s -> entrance pupil diam (next shell)= %s"%(d,dNext)
        dNext=dNext- (1. *2)
        print "inner (without walls)= %s"%(dNext)
        print "\tintersection diam (next shell)= %s"%inverso(F,dNext,H[1])
        print 10*"-"


F=1000.
D=346.807377
SH=[0,300.]
tpars=[0.0008578 ,0]
#print diretto(F,D,300.,wolter=1)
for i in range(0,12): print i,Nouter(i,D,F,SH,tpars,ffdeg=0.1)
        