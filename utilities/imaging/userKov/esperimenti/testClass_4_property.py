# schema proprieta' get e set ##############
class C(object):
    def __init__(self): self.__x = None
    def getx(self): return self.__x
    def setx(self, value): self.__x = value
    def delx(self): del self.__x
    x = property(getx, setx, delx, "I'm the 'x' property.")
# ###########################################

class cc(object):
    def __init__(self,R):
        self.R = R
    def getx(self):
        print "getx ",self.R
        return self.R*2
    def setx(self, value):
        print "setx ",value
        self.R = value/2
    D = property(getx, setx, "I'm the 'D' property.")

class cc2:
    def __init__(self,R):
        self.R = R
    def getx(self):
        print "getx ",self.R
        return self.R*2
    def setx(self, value):
        print "setx ",value
        self.R = value/2
    D = property(getx, setx, "I'm the 'D' property.")

if __name__=="__main__":
    a=cc(4)
    b=cc2(4)
    print "\n\na=cc(4)  eredita da object"
    print "b=cc2(4)  uguale ad a, ma non eredita\n"
    print "\n---------------"    
    print "a.R: ",a.R
    print "a.D:\n",a.D,"\n"
    print "b.R: ",b.R
    print "b.D:\n",b.D
    print "---------------\n"
    print "assegnando nuovi valori:"
    print "a.D=6"
    a.D=6    
    print "a.R= ",a.R,"\n"
    print "b.D=6 (setx non viene chiamata..)"
    b.D=6
    print "b.R= ",b.R
    


    