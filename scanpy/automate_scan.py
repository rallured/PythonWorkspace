#Methods to automate data acquisition for testing the proportional counter on
#the BRP project in the Kaaret group at U of Iowa.

#Author: Sean Mattingly
#Date:   Jan 27th, 2011

import mca_mod
import thorALL
import os
import time
from scipy import *
from numpy import *


def test_range_auto(x_bounds, y_bounds, x_jog_dist = 1, y_jog_dist = 1,
                    ac_time = 60):
    """Function to automate acquiring data over a rectangle as given by the
    manually - measured left, right, top and bottom constraining y and x
    values.
    
    All instances of left and right are if one is looking down the tube in the
    same directions the xrays come out.
    
    As written, this program requires one to give the ABSOLUTE positions of
    the desired square in terms of the THOR mount coordinates. E.G. this
    program will NOT get the current position and do that work for you!

    Another note: It is possible that the THOR mount may have the x or y
    coordinate inverted - take this into account when giving the bounds on the
    desired scanning bounds. Generally, if one just manually finds the values
    and throws those in from the THOR control panel, those will work.
    
    Required inputs:
    x_bounds:
        Tuple (x1, x2) of Floats of length 2. Bounds delimiting the x range to 
        scan over.
        Order doesn't matter.
    y_bounds:
        Tuple (y1, y2) of floats of length 2. Bounds delimiting the y range to 
        scan over.
        Order doesn't matter.
        
    Optional Inputs:
    x_jog_dist:
        Int. Value to jog the x coordinate by while scanning.
    y_jog_dist:
        Int. Value to jog the y coordinate by while scanning.
    ac_time:
        Int. Number of seconds to acquire data for.
        
    
    """
    folder_name = time.asctime(time.localtime()).replace(' ', '_').replace(':','_')
    #Platform independence using os
    os.mkdir(folder_name)
    #Initialize THOR mount connection.
    thorALL.init()
    raw_input('Change the settings if it needs to have a min pos < 0 or a max'
               + ' pos > 100, then press enter: ')
    #Initialize MCA
    mca_mod.init()
    #Scan over it all!
    for pos in gen_arrays(x_bounds, y_bounds, x_jog_dist, y_jog_dist):
        #Move to the current position
        print('Moving to ' + str(pos))
        move(pos)
        #Wait to attempt to alleviate strange data issue.
        time.sleep(20)
        #Take data
        data = mca_mod.acquire(ac_time,disc=10)
        #Write the data
        fname_root = folder_name + os.sep + 'spec'
        suf = 'x_' + str(pos[0]) + 'y_' + str(pos[1])
        mca_mod.save(fname_root+suf+'.fits',keywords=[('XPOS',pos[0],'X position'),('YPOS',pos[1],'Y position')])
    thorALL.close()
    return

def gen_arrays(x_bounds, y_bounds, x_jog_dist, y_jog_dist):
    """Generates a list of (x,y) coordinates to iterate through
    """
    coords = []
    #Obviously, this isn't the most optimal way  
    #to organize these loops, but it
    #should be intuitive and readable.
    for x in arange(min(x_bounds), max(x_bounds) + x_jog_dist, step=x_jog_dist):
        for y in arange(min(y_bounds), max(y_bounds) + y_jog_dist, step=y_jog_dist):
            coords.append((x, y))
    return coords
            
def move(pos):
    """Moves the THOR mount to a given x, y position given by the tuple pos = (x,y).
    
    Required inputs:
    pos:
        Float Tuple of length 2. Corresponds to the (x, y) position of the desired position for the THOR mount.
    """
    thorALL.h_move(pos[0])
    thorALL.v_move(pos[1])
    return
    
    



























