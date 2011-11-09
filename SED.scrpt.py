#!/usr/bin/env python
# 
#  SED.scrpt.py
#  Simulation Software
#  
#  Created by Alexander Rudy on 2011-11-09.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 



import math, copy, sys, time, logging, os, argparse

LongHelp = """
This is the Command Line Interface to the SEDMachine Simulator Program.

The program works in the following stages:
    1. Instrument Property Setup
    2. Spectrum Generation
    3. Spectra Sub-Image Generation
    4. Full Image Generation (Placing Sub-Images)
    5. Noise Image Generation
    6. Scattered Light Generation
"""

parser = argparse.ArgumentParser(description="Utilities for simulating the SEDMachine",epilog=LongHelp,formatter_class=argparse.RawDescriptionHelpFormatter,)
parser.add_argument('--no-cache',action='store_const',const=True,default=False,help="Re-run all algorithm parts")
parser.add_argument('--stage',action='store',metavar='N',nargs=1,type=int,help="Start at stage N")
parser.add_argument('-t','--thread',action='store_const',const=True,default=False,help="Enable Threading")
parser.add_argument('-n',action='store',metavar='N',nargs=1,type=int,help="Limit to the first N spectra")
parser.add_argument('-o',action='store',metavar='N',nargs=1,type=int,help="Offset the starting spectra number by N")
parser.add_argument('-d',action='store_const',const=True,default=False,help="Enable Debugging Mode")
specgroup = parser.add_argument_group('Spectra')
specgroup.add_argument('-s',choices='bf',help="Select Spectrum: b: BlackBody f:Flat")
specgroup.add_argument('-T',action='store',metavar='Temp',type=float,default=5000.0,help="BlackBody Temperature to use")
specgroup.add_argument('-V',action='store',metavar='Value',type=float,default=1000.0,help="Flat Spectrum Value")
args = parser.parse_args()

class Simulation(object):
    """This is the master class for a complete simulation"""
    def __init__(self):
        super(Simulation, self).__init__()
        self.stages = {}
        
        
    def add_stage(self,stage):
        """docstring for add_stage"""
        label = stage.label
        self.stages[label] = stage
        
    def run_simulation(self,starting_stage=1):
        """docstring for run_simulation"""
        NStages = len(self.stages)
        AStages = range(starting_stage,NStages)
        print AStages
        
        
        
        
        
        
        
        