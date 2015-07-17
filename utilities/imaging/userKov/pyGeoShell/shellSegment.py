class shellSegment(object):
    '''segment of shell with cylindrical symmetry.
        input arguments:
            length: length of shell projection along z
            zOffset:
            angularSize:
            profile|polypars|slope'''

    def __init__(self,length,zOffset=0.,angularSize=2.,profile=None,polyCoeff=None,slope=None):
        self.length=length
        if polyCoeff is not None:
            self.polyCoeff=polyCoeff
            self.calculateProfile

    def calculate_profile(self,z):
        '''calculate the profile at a given z from the polynomial coefficients
        sorted from zeroth grade.'''
        