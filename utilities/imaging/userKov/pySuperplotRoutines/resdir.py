class resDir :
    def __init__(self,ar,ndim=0):
        self.data=array(ar)
        self.ndim=ndim
        if ndim==0:
            self.ndim=len (ar[0])-1
            print "numero di dimensioni non dato!"
            print "fissato in" , self.ndim
            #print dir(self)
    def __getitem__ (self,index):
        return self.data[index]
    def __repr__(self):
        s=""
        for res in self.data:
            for el in res:
                s=s+(str(el)+"\t")
            s=s+("\n")
        return s
    def sel(self,nd=-1):
        if nd==-1:nd=self.ndim+1
        bd=self.data[0::nd]
        n=resDir(bd,nd)
        return n
    def put(self,file):
        f=open(file,"w")
        f.write(str(self))
        f.close()