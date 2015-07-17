'''esperimenti di subclassamento e wrapper,
generando output nella routine __getattr__.
permette di vedere che metodi speciali vengono chiamati.
E' interessante notare che viene richiamata anche
quando nel prompt si attiva il completamento automatico dopo il .
in questo caso chiama: (__members__ e  __method__)
si possono vedere altri casi (trucchi per non far casino con il
completamento automatico: copia e incolla, history o scrivere lettera
in piu' poi cancellarla):

print a6 ->  __str__
a6 (nel prompt) -> __repr__
'''

class resList5(list):
    '''eredita da lista'''
    def __init__(self,data=[]):
        list.__init__(self,data)
        
    def __getattr__(self,attr):
        print  "\nin __getattr__: attr= ",attr
        #b= getattr(self,attr) se usavo questo ottenevo ricorsione infinita
        #perche' richiamava metodo di resList5, non di list
        list.__getattr__(self,attr)
        print "\nb= getattr(self.data,attr),"
        print "print b,self.data:",b,self
        print "return b:"
        return b

    def prova(self):
        print "scrivo 'prova'!!"
        
class resList6:
    ''''voglio provare a ottenere un effetto simile senza ereditarieta'.
    per gli array infatti non si puo' ereditare'''
    def __init__(self,data=[]):
        self.data=list(data)
        
    def __getattr__(self,attr):
        print  "\nin __getattr__: attr= ",attr
        print "print self.data:",self.data
        b= getattr(self.data,attr)
        print "b= getattr(self.data,attr)"
        print "print b,self.data:",b,self.data
        print "return b:"
        return b

    def prova(self):
        print "scrivo 'prova'!!"

if __name__=="__main__":
    print "\n-------------------------------------"*2
    print "a5=resList5([1,-2,3]) [subclassa lista]"
    print "a6=resList6([1,-2,3]) [non subclassa]"
    print "---------------------------------------"
    print "---------------------------------------"
    a6=resList6([1,-2,3])
    a5=resList5([1,-2,3])
    print "\n---- a5, subclassato ------"
    print "\n- a5: ",a5
    print "\n- a5.sort: ",a5.sort   #il metodo esiste (subclassato)
    print "\n- a5: ",a5
    print "\n- a5.sort(): ",a5.sort()    #ma non restituisce valore perche' non viene richiamato espressamente
    print "\n- a5: ",a5
    try:
        print "\n- a5.kov (proprieta' inesistente):",a5.kov
    except:
        print "riscontra errore\n"
    print "\n- a5.prova(): ",a5.prova()
    print "---- a6, non subclassato --"
    print "- a6.data: ",a6.data
    print "\n- a6.sort: ",a6.sort
    print "\n- a6.data: ",a6.data
    print "\n- a6.sort(): ",a6.sort()
    print "\n- a6.data: ",a6.data
    try:
        print "\n- a6.kov (proprieta' inesistente):",a6.kov
    except:
        print "riscontra errore"
    print "\n- a6.prova(): ",a6.prova()