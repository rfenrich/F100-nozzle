# -*- coding: utf-8 -*-
"""
Material module for axisymmetric nozzle

Rick Fenrich 6/28/16
"""

class Material:
    def __init__(self, k, alpha, E, v):
        self.k = k # W/m*K, thermal conductivity of wall
        self.alpha = alpha # 1/K, coeff. of thermal expansion
        self.E = E # Pa, elastic modulus
        self.v = v # Poisson's ratio
    
