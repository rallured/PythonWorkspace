'''comportamento delle proprieta' assegnate per riferimento'''
from numpy import *
 
class oggetto:
    '''
    '''
    def __init__(self,data):
        self.ar=array(data)
    def add1(self):
        print "la funzione add1 modifica a.ar inserendo in a.ar valori in sequenza"
        print "l'effetto sparisce all'uscita dalla routine."
        for i in self.ar:
            i=i+1
    def add2(self):
        print "la funzione add2 modifica a.ar creando un nuovo array"
        print "con elementi in sequenza ponendolo nella proprieta' self.ar.\n"
        copia=[]
        for i in self.ar:
            copia.append(i+1)
        self.ar=array(copia)


if __name__=="__main__":
    a=oggetto([1,-2,4])
    b=array([4,-1,2])
    a.ar=b
    print "\n\n-------------------------"
    print "a.ar= %s | b= %s"%(a.ar,b)
    print "\n-eseguo a.add1"
    a.add1()
    print "a.ar= %s | b= %s"%(a.ar,b)
    print "\n- eseguo a.add2"
    a.add2()
    print "a.ar= %s | b= %s"%(a.ar,b)