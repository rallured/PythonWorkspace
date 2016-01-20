import numpy as np
import matplotlib.pyplot as plt

def inducedR(stress,young,poisson,hs,hf):
    """Return the induced radius of curvature on a plate-like
    substrate with a given Young's modulus, Poisson ratio,
    substrate thickness, film thickness, and coating stress
    stress and young in units of MPa
    hf,hs, and returned value in same length unit
    """
    young = young/(1-poisson)
    return -young*hs**2/6/hf/stress

def sag(diam,R):
    """Return the peak-to-valley sag from a bow over a given
    diameter, both in same units"""
    return diam**2/8./R

def inducedBow(stress,young,poisson,hs,hf,diam):
    """Return the induced peak-to-valley bow induced on a
    circular plate. Calls inducedR and then uses the sag
    equation.
    """
    R = inducedR(stress,young,poisson,hs,hf)
    return sag(diam,R)
