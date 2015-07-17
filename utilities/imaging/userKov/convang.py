import math


def r2d(an):return an/(math.pi/180.)
def d2r(an):return an*(math.pi/180.)
def s2r(an):return d2r(an)/3600.
def r2s(an):return r2d(an)*3600.
def m2r(an):return d2r(an)/60.
def r2m(an):return r2d(an)*60.

if __name__=='__main__':
    print 'conversion function for angles:'
    print 'r2d(an) radians to degrees'
    print 'r2s(an) radians to arcsec'
    print 'r2m(an) radians to arcmin'
    print 'and viceversa\n'
    print 'example: r2m(2.0) ->', r2m(2.0)