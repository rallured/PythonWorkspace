'''1- confronto tra oggetti
2- descrittori
'''

###1-confronto tra classi
class P(object):
    def __init__(self,n):
        self.nome=n

class Persona(object):
    def __init__(self,n):
        self.nome=n
    def __eq__(self,other):
        return self.nome == other.nome

###2-descrittori
class desc(object):
    def __init__ (self,initval=None,name='var'):
        self.val=initval
        self.name=name
    def __get__(self,obj,objtype):
        print "get: ",self.name
        return self.val
        print "obj: ",obj
        print "self: ",self
    def __set__(self,obj,val):
        print "set: ",self.name
        print self.val
        print "obj: ",obj
        print "self: ",self

class my(object):
    x=desc(10,'var "x"')
    y=5    

if __name__=="__main__":
###1-confronto tra classi
    #con una classe che comprende il metodo __eq__:
    kov=Persona("Vin")
    kov2=Persona("Vin")
    kov3=Persona("Vince")
    kov4=kov
    print "####################"
    print "kov: ",kov
    print "kov.nome= ",kov.nome    
    print "kov2: ",kov2
    print "kov2.nome= ",kov2.nome    
    print "kov3: ",kov3
    print "kov3.nome= ",kov3.nome
    print "kov4: ",kov4
    print "kov4.nome= ",kov4.nome    
    print "kov==kov2, kov==kov3, kov==kov4:\n",kov==kov2, kov==kov3, kov==kov4
    print "####################"    
    #con una classe che non lo comprende:
    kov=P("Vin")
    kov2=P("Vin")
    kov3=P("Vince")
    kov4=kov
    print "\n\n####################"
    print "kov: ",kov
    print "kov.nome= ",kov.nome    
    print "kov2: ",kov2
    print "kov2.nome= ",kov2.nome    
    print "kov3: ",kov3
    print "kov3.nome= ",kov3.nome
    print "kov4: ",kov4
    print "kov4.nome= ",kov4.nome    
    print "kov==kov2, kov==kov3, kov==kov4:\n",kov==kov2, kov==kov3, kov==kov4
    print "####################"
    
###2-descrittori    
