# -*- coding: utf-8 -*-
"""
Material module for axisymmetric nozzle

Rick Fenrich 6/28/16
"""

class Inlet:
    def __init__(self, PstagInlet, TstagInlet):
        self.Pstag = PstagInlet # Pa, stagnation pressure of inlet
        self.Tstag = TstagInlet # K, stagnation temp. of inlet
        
    def setMach(self, mach):
        self.mach = mach