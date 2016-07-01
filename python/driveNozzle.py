# -*- coding: utf-8 -*-
"""
Test driver script for nozzleIdeal.py and nozzleNonIdeal.py

Provides inputs to the nozzleIdeal and nozzleNonIdeal functions, runs them, and
produces plots comparing ideal and non-ideal nozzle behavior, as well as non-
ideal nozzle heating and friction characteristics. Based on the equivalent
Matlab code written from 7/29/15 onwards.

Rick Fenrich 6/27/16
"""

import numpy as np
import nozzle
import geometry
import material
import component
import inlet
import environment
import fluid
import mission

import quasi1D

config = {}; 

# ========================== INPUT PARAMETERS ================================
config["mission"] = 3 # mission number
config["fidelity"] = "med" # "low" (quasi-1D flow), "med" (Euler), or "high" (RANS)

if(config["fidelity"] == "med"):
    config["meshSize"] = "coarse" # "coarse", "medium", or "fine"
    config["governing"] = "euler" # "euler" or "rans"
    
# ------------------------ SET HEAT TRANSFER PARAMS --------------------------
config["hInf"] = 25 # W/m^2/K, heat transfer coeff. from external wall to env.

# ------------------------ SET MATERIAL PROPERTIES ---------------------------
config["thermalConductivity"] = 8.6 # W/m*K, thermal conductivity of wall
config["coeffThermalExpansion"] = 2.3e-6 # 1/K, coeff. of thermal expansion
config["elasticModulus"] = 80e9 # Pa, elastic modulus
config["poissonRatio"] = 0.3 # Poisson's ratio

# --------------------------- SET FLIGHT REGIME ------------------------------
if(config["mission"] == 1): # static sea-level thrust case
    altitude = 0.
    mach = 0.01
    config["inletStagTemp"] = 888.3658
    config["inletStagPres"] = 3.0550e5
elif(config["mission"] == 2): # intermediate case
    altitude = 15000.
    mach = 0.5
    config["inletStagTemp"] = 942.9857
    config["inletStagPres"] = 2.3227e5
elif(config["mission"] == 3): # high speed, high altitude case
    altitude = 35000.
    mach = 0.9
    config["inletStagTemp"] = 1021.5
    config["inletStagPres"] = 1.44925e5
elif(config["mission"] == 4): # case with shock in nozzle
    altitude = 0.
    mach = 0.01
    config["inletStagTemp"] = 900.
    config["inletStagPres"] = 1.3e5
elif(config["mission"] == 5): # subsonic flow
    altitude = 0.
    mach = 0.01
    config["inletStagTemp"] = 900.
    config["inletStagPres"] = 1.1e5

# --------------------- SET ERROR TOLERANCE RANGES ---------------------------
# Set err tolerances for various iterations and solvers
if(config["fidelity"] == "low"):
    config["errorExitTempBetweenIterations"] = 1e-8
    config["errorSolverApparentThroatLocation"] = 1e-6
    config["errorSolverRelativeMachSquared"] = 1e-10
    config["errorSolverAbsoluteMachSquared"] = 1e-10
    config["errorDenominatordMdx"] = 4; # this is not an error tolerance, 
        # rather it is used to set the slope of dMdx in the transonic regime

# ------------------------ SET NOZZLE GEOMETRY -------------------------------
config["innerWallParameterization"] = "B-spline" # options include 'linear', 'spline', 'B-spline', and 'B-spline-mex'

# SKIP IMPLEMENTATION OF CUBIC SPLINE GEOMETRY

if(config["innerWallParameterization"] == "B-spline"):
    config["bSplineCoefficients"] = np.array(([0.0000, 0.0000, 0.1500, 0.1700, 
        0.1900, 0.2124, 0.2269, 0.2734, 0.3218, 0.3218, 0.3230, 0.3343, 0.3474, 
        0.4392, 0.4828, 0.5673, 0.6700, 0.6700],[0.3255, 0.3255, 0.3255, 
        0.3255, 0.3255, 0.3238, 0.2981, 0.2817, 0.2787, 0.2787, 0.2787, 0.2797, 
        0.2807, 0.2936, 0.2978, 0.3049, 0.3048, 0.3048]))
        
# ---------------------- SET NOZZLE WALL GEOMETRY ----------------------------
config["thicknessParameterization"] = "piecewise-linear"

if(config["thicknessParameterization"] == "piecewise-linear"):
    config["thicknessNodeArray"] = np.array(([0., 0.33, 0.67],[0.01, 0.01,   \
      0.01]))
      
# -------------------- SET NOZZLE FLUID PROPERTIES ---------------------------
config["specificHeatRatio"] = 1.4
config["specificGasConstant"] = 287.06

# ====================== SET UP NOZZLE DEFINITION ============================
nozzle = nozzle.Nozzle()

nozzle.wall = component.AxisymmetricWall()

nozzle.wall.geometry = geometry.Bspline(config["bSplineCoefficients"])
#nozzle.wall.geometry = geometry.PiecewiseLinear(np.array(([0., 0.67],[0.32, 0.24])))
nozzle.wall.thickness = geometry.PiecewiseLinear(config["thicknessNodeArray"])
nozzle.wall.material = material.Material(config["thermalConductivity"],      \
  config["coeffThermalExpansion"], config["elasticModulus"],                 \
  config["poissonRatio"])

#nozzle.baffles = VerticalBaffles()
#nozzle.stringers = AxisymmetricStringers()

nozzle.inlet = inlet.Inlet(config["inletStagPres"],config["inletStagTemp"])
#nozzle.inlet.setMach(0.2)

nozzle.environment = environment.Environment(altitude,config["hInf"])

nozzle.fluid = fluid.Fluid(config["specificHeatRatio"],                      \
  config["specificGasConstant"])
  
nozzle.mission = mission.Mission(config["mission"])
nozzle.mission.setMach(mach)

# ========================== RUN CALCULATIONS ================================

tol = {}
tol["exitTempPercentError"] = 1e-8
tol["solverRelTol"] = 1e-10
tol["solverAbsTol"] = 1e-10
tol["solverApparentThroatLocation"] = 1e-6
results = quasi1D.analysis(nozzle,tol)
  


        
    
    