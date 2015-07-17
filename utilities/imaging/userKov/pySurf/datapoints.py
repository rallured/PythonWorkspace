import numpy as np

class DataPoints(np.ndarray):
    """a set of datapoints."""
    def __init__(self):
        pass
        
    def read(self, filename,**kwargs):
        self.points=np.genfromtxt(filename,**kwargs)
    
    def save(self,filename,binary=false):
        np.savetxt(filename, X, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ')
    
    

if __name__=="__main__":
    filename=
    markersFile=
    pp=np.genfromtxt(filename)
    m1,m2=readmarkers(markersFile)
    trans=find_affine(m1,m2)
    trans(points)
    save(points)
    

