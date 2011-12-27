# 
#  Utilities.py
#  Simulation Software
#  
#  Created by Alexander Rudy on 2011-11-02.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 

import numpy as np
import pyfits as pf
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt


def bin(array,factor):
    """Bins an array by the given factor"""
    
    finalShape = tuple((np.array(array.shape) / factor).astype(np.int))
    Aout = np.zeros(finalShape)
    
    for i in range(factor):
        Ai = array[i::factor,i::factor]
        Aout += Ai
    
    return Aout


def expandLim(axis,scale=0.05):
    """Expands Axis Limits by *scale*"""
    xmin,xmax,ymin,ymax = axis
    xran = abs(xmax-xmin)
    yran = abs(ymax-ymin)

    xmax += xran*scale
    xmin -= xran*scale
    ymax += yran*scale
    ymin -= yran*scale

    axis = (xmin,xmax,ymin,ymax)
    return axis

def rangemsg(array,name):
    """Message describing this array"""
    MSG = "Array named %(name)s has %(elements)d els with shape %(shape)s. Range %(range)s. Zeros %(zeros)d (%(zper)3d%%). NaNs %(nans)d (%(nper)3d%%). Type %(type)3s"
    fmtr = {}
    fmtr["elements"] = array.size
    fmtr["shape"] = str(array.shape)
    fmtr["name"] = name
    fmtr["min"] = np.min(array)
    fmtr["max"] = np.max(array)
    fmtr["zeros"] = np.sum(array == np.zeros(array.shape))
    fmtr["zper"] = float(fmtr["zeros"]) / float(fmtr["elements"]) * 100
    fmtr["nans"] = np.sum(np.isnan(array))
    fmtr["nper"] = float(fmtr["nans"]) / float(fmtr["elements"]) * 100
    fmtr["type"] = array.dtype
    fmtr["range"] = "[%(min)5.5g,%(max)5.5g]" % fmtr
    return MSG % fmtr
