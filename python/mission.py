# -*- coding: utf-8 -*-
"""
Mission module for axisymmetric nozzle

Rick Fenrich 6/28/16
"""

class Mission:
    def __init__(self, number):
        self.number = number # number of mission
        
    def setMach(self, mach):
        self.mach = mach # mach number
    
