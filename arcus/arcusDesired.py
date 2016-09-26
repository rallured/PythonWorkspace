import numpy as np

def desiredGratRef(indexYaw,indexRoll,indexPitch,\
                   shiftRoll,shiftPitch,shiftYaw,gratNum):
    """index is from theodolite (only pitch from far mirror)
    theodolite measurements are from reference grating
    shift is from WFS and laser and is from pre-bond previous grating
    Returns yaws in arcsec and pitch/roll in WFS units
    gratNum is number of gratings from the reference grating"""
    yaw = indexYaw - shiftYaw
    roll = -indexRoll - shiftRoll/.034
    pitch = (51.4*gratNum)-indexPitch - shiftPitch/.034

    print 'Referenced to Previous Grating:'
    print 'Grat Yaw: ' + str(yaw)
    print 'Grat Roll: ' + str(roll*.034)
    print 'Grat Pitch: ' + str(pitch*.034)

    yaw = indexYaw
    roll = -indexRoll
    pitch = (51.4*gratNum)-indexPitch

    print 'Absolute Reference:'
    print 'Abs Yaw: ' + str(yaw)
    print 'Abs Roll: ' + str(roll*.034)
    print 'Abs Pitch: ' + str(pitch*.034)
    
    return

def desiredAbsRef(theoYaw,indexRoll,indexPitch):
    """theo is from theodolite, index is from WFS,
    shift is from WFS and laser"""
    yaw = theoYaw
    roll = indexRoll/.034
    pitch = indexPitch/.034
    return roll, pitch, yaw
