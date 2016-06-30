# -*- coding: utf-8 -*-
"""
Fluid module

Rick Fenrich 6/28/16
"""

import numpy as np

class Fluid:

    Tinterp = np.array(([175, 200, 225, 250, 275, 300, 325, 350, 375, 400,   \
      450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050,     \
      1100, 1150, 1200, 1250, 1300, 1350, 1400, 1500]))    
    
    def __init__(self, gamma, R):
        self.type = "air"
        self.gam = gamma # ratio of specific heats
        self.R = R # J/kg/K specific gas constant
        
    def Pr(self, T): # Prandtl number
        PrNum = np.array(([0.744, 0.736, 0.728, 0.72, 0.713, 0.707, 0.701,   \
          0.697, 0.692, 0.688, 0.684, 0.68, 0.68, 0.68, 0.682, 0.684, 0.687, \
          0.69, 0.693, 0.696, 0.699, 0.702, 0.704, 0.707, 0.709, 0.711,      \
          0.713, 0.715, 0.717, 0.719, 0.722]))
        return np.interp(T,self.Tinterp,PrNum)
        
    def k(self, T): # thermal conductivity
        thermalConductivity = 0.01*np.array(([1.593, 1.809, 2.02, 2.227,     \
          2.428, 2.624, 2.816, 3.003, 3.186, 3.365, 3.71, 4.041, 4.357,      \
          4.661, 4.954, 5.236, 5.509, 5.774, 6.03, 6.276, 6.52, 6.754, 6.985,\
          7.209, 7.427, 7.64, 7.849, 8.054, 8.253, 8.45, 8.831]))
        return np.interp(T,self.Tinterp,thermalConductivity)
        
    def Cp(self,T): # specific heat at constant pressure
        specificHeatCP = np.array(([1002.3, 1002.5, 1002.7, 1003.1, 1003.8,  \
          1004.9, 1006.3, 1008.2, 1010.6, 1013.5, 1020.6, 1029.5, 1039.8,    \
          1051.1, 1062.9, 1075.0, 1087.0, 1098.7, 1110.1, 1120.9, 1131.3,    \
          1141.1, 1150.2, 1158.9, 1167.0, 1174.6, 1181.7, 1188.4, 1194.6,    \
          1200.5, 1211.2]))
        return np.interp(T,self.Tinterp,specificHeatCP)
    
