import numpy as np
import matplotlib.pyplot as plt

def sphericalCalc(colToCGH=100.,colF=400.,fieldF=200.):
    #Compute location of spherical image after collimator
    colImg = 1./(1./colF - 1/(colToCGH+420.))
    colMag = colImg/(colToCGH+420.)
    print 'Image at %.2f after collimator with %.2f mag' % (colImg,colMag)
    #Compute required location of field lens
    mag = 14.6/25.4/colMag
    print 'Field lens needs %.2f magnification' % mag
    fieldObj = fieldF/mag-fieldF
    fieldToCol = colImg-fieldObj
    print 'Field lens at %.2f after collimator' % fieldToCol
    fieldImg = 1./(1/fieldF+1/fieldObj)
    print 'Final image is %.2f after field lens' % fieldImg

    return fieldToCol,fieldImg

def cylCalc(fieldToCol,cylF=1000.,colToCGH=100.,colF=400.,fieldF=200.):
    #Compute location of image after CGH
    cghImg = 1./(1./200.-1./420.)
    cghMag = cghImg/420.
    colObj = colToCGH-cghImg
    colImg = 1./(1./colF-1./colObj)
    colMag = colImg/colObj * cghMag
    print 'Image at %.2f after collimator with %.2f mag' % (colImg,colMag)
    fieldObj = fieldToCol-colImg
    fieldImg = 1./(1./fieldF-1./fieldObj)
    fieldMag = fieldImg/fieldObj * colMag
    print 'Image at %.2f after field lens with %.2f mag' % (fieldImg,fieldMag)
    #Determine where cylindrical field lens should go
    
