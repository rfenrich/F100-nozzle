# -*- coding: utf-8 -*-
"""
Fluid module

Rick Fenrich 6/28/16
"""

class Fluid:
    def __init__(self, gamma, R):
        self.type = "air"
        self.gam = gamma # ratio of specific heats
        self.R = R # J/kg/K specific gas constant
    
