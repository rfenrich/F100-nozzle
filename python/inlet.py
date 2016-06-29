# -*- coding: utf-8 -*-
"""
Material module for axisymmetric nozzle

Rick Fenrich 6/28/16
"""

class Inlet:
    def __init__(self, PstagInlet, TstagInlet):
        self.PstagInlet = PstagInlet, # Pa, stagnation pressure of inlet
        self.TstagInlet = TstagInlet # K, stagnation temp. of inlet