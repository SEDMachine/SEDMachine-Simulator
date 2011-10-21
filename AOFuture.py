# 
#  AOFuture.py
#  Simulation Software
#  
#  Created by Alexander Rudy on 2011-10-19.
#  Copyright 2011 Alexander Rudy. All rights reserved.
#  File contains functions which will be included in the next release of AstroObjects

from AstroObject import *

import logging


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


def disable_Console():
    """Disables console Logging"""
    logging.getLogger('').removeHandler(console)

def enable_Console():
    """docstring for enable_Console"""
    logging.getLogger('').addHandler(console)