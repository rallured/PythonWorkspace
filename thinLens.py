import numpy as np
import matplotlib.pyplot as plt

def computeImg(F,obj):
    """
    Compute the image location and magnification
    based on lens focal length and object position.
    """
    img = 1./(1./F-1./obj)
    mag = img/obj
    return img,mag

