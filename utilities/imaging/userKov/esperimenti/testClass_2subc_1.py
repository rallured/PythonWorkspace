'''tentativi di subclassare array:
-kovint: subclassamento di intero, funzionante.
-resList: subclassamento di lista, funzionante. '''

###1-subclassamento di intero
class kovint(int):
    #subclassamento di intero funzionante, con tutti metodi e proprieta'
    #in piu' ha la proprieta' val che contiene il valore e che e' l'unico modo
    #di modificarlo 
    def __init__(self,v):
        self.val=v
        int.__init__(self)

if __name__=="__main__":
    print "\n######################################"
    print "#### Subclassamento di intero: #########\n"
    a=kovint(5)
    print "a=kovint(5)\na, a.val: ",a,a.val
    a.val=7
    print "\na.val=7"
    print "a, a.val: ",a,a.val
    a=8
    print "\na=8"
    try:
        print "a, a.val: ",a,a.val   
    except:
        #print "a, a.val: ",a
        print "\na e' un int, non ha attributo val!"
        print "\n###################################\n\n"

###2-subclassamento di lista
from numpy import *

class resList(list):
    #questa funziona (credo), subclassa la lista e ne eredita i metodi
    #a differenza del precedente non puo' essere modificato, perche' __data
    #e' di sola lettura per via dei trattini iniziali
    def __init__(self,data=[]):
        #print data
        self.__data=data
        list.__init__(self,self.__data)

class resList2(list):
    # anche questa funzionerebbe, in modo ancora piu' semplice,
    # per accedere ai dati nel codice della classe si accede a self
    # che si comporta come una lista
    def __init__(self,data=[]):
        list.__init__(self,data)
    def esempio(self):
        print self[2]
        
if __name__=="__main__":
    print "\n######################################"
    print "#### Subclassamento di lista: #########\n"
    a=resList([1,2,3])
    print"\na=resList([1,2,3]) \na, id(a), type(a): ",a, id(a), type(a)
    a=resList([4,5,6])
    print"\na=resList([4,5,6]) \na, id(a), type(a): ",a, id(a), type(a)
    a=[4,5,6]
    print"\na=[4,5,6] \na, id(a), type(a): ",a, id(a), type(a)
    a=resList([4,5,6])
    try:
        print "a=resList([4,5,6]) \na.__data: ",a.__data
    except:
        print "_\n_data e' una varibile privata,\n python da' errore."
    b=resList2([2,4,6])
    print "\n\nb=resList2([2,4,6])"
    print "b.esempio() -> stampa il terzo valore:\n"
    b.esempio()
    print "\n###############################"
 