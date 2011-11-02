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