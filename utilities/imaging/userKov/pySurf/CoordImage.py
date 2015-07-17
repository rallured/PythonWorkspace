import numpy as np


class CoordImage(np.ndarray):
    
    def __init__(self, xgrid=None, ygrid=None, *args, **kwargs):
        if xgrid!=None:self.xgrid=xgrid
        if ygrid!=None:self.ygrid=ygrid
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self
    
    def __setattr__(self, k, v):
        self.__dict__[k] = v
    def __getattr__(self, k):
        ''' we don't need a special call to super here because getattr is only 
            called when an attribute is not found in the instance's dictionary'''
        return self.data[k]
    
