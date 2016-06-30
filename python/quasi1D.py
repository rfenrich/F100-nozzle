# -*- coding: utf-8 -*-
"""
Perform quasi-1D area-averaged Navier-Stokes analysis on axisymmetric nozzle

Rick Fenrich 6/28/16
"""

import numpy as np
import scipy.optimize
import scipy.integrate
import time

#==============================================================================
# Sutherland's Law of dynamic viscosity of air
#==============================================================================
def dynamicViscosity(T):
    mu = 1.716e-5*(T/273.15)**1.5*(273.15 + 110.4)/(T + 110.4)
   
#==============================================================================
# Area-Mach function from 1-D mass conservation equations
#==============================================================================
def areaMachFunc(g,M):
    ((g+1)/2)**((g+1)/(2*(g-1)))*M/(1+(g-1)*M**2/2)**((g+1)/(2*(g-1)))

#==============================================================================
# Ideal analysis of nozzle (no heat transfer or friction)
#==============================================================================
def analysisIdeal(nozzle):      
    pass

#==============================================================================
# Determine state of nozzle stagnation pressure and temperatures at throat and
# exit, nozzle geometry, and pressure ratio between reservoir and atmosphere
#==============================================================================
def nozzleState(nozzle,pressureRatio,PsT,TsT,PsE,TsE):
    gam = nozzle.fluid.gam
    Athroat = nozzle.wall.geometry.area(nozzle.wall.geometry.xThroat)
    Aexit = nozzle.wall.geometry.area(nozzle.wall.geometry.length)
    
    rootFunc = lambda x: ((gam+1)/2)**((gam+1)/(2*(gam-1)))*x/(1+(gam-1)*    \
      x**2/2)**((gam+1)/(2*(gam-1)))*np.sqrt(TsT/TsE)*(PsE/PsT) - Athroat/Aexit
      
    MsubsonicCritical = scipy.optimize.fsolve(func=rootFunc,x0=0.5,xtol=1e-12)
    MsupersonicCritical = scipy.optimize.fsolve(func=rootFunc,x0=2.,xtol=1e-12)
    
    PtRatioSubsonic = (1 + (gam-1)*MsubsonicCritical**2/2)**(gam/(gam-1))
    PtRatioSupersonic = (1 + (gam-1)*MsupersonicCritical**2/2)**(gam/(gam-1))
      
    MbehindShock = np.sqrt((1 + (gam-1)*MsupersonicCritical**2/2)/           \
      (gam*MsupersonicCritical**2 - (gam-1)/2))
    normalShockPtRatio = (Aexit/Athroat)*((gam+1)/2)**((gam+1)/(2*(gam-1)))  \
      *MbehindShock*np.sqrt(1 + (gam-1)*MbehindShock**2/2)
      
    deltaPtRatio = 0.05 # a tolerance to determine if fully expanded flow occur
    Mshock = 0.
    if( pressureRatio <= 1 ):
        status = "no flow"
    elif( pressureRatio < PtRatioSubsonic ):
        status = "subsonic"
        Mexit = np.sqrt(2/(gam-1))*np.sqrt(pressureRatio**((gam-1)/gam) - 1)
        Mshock = Mexit # return Mexit in place of Mshock
    elif( pressureRatio < normalShockPtRatio):
        status = "shock"
        
        rootFunc2 = lambda x: ((gam+1)/2)**((gam+1)/(2*(gam-1)))*x*          \
          np.sqrt(1 + (gam-1)*x**2/2) - pressureRatio*(Athroat/Aexit)*       \
          np.sqrt(TsE/TsT)
        Mexit = scipy.optimize.fsolve(func=rootFunc2,x0=0.5,xtol=1e-12)
        PtRatio = (Athroat/Aexit)*np.sqrt(TsE/TsT)/(((gam+1)/2)**((gam+1)/   \
          (2*(gam-1)))*Mexit/(1+(gam-1)*Mexit**2/2)**((gam+1)/(2*(gam-1))))
        
        rootFunc3 = lambda x: (((gam+1)*x**2/2)/(1 + (gam-1)*x**2/2))**      \
          (gam/(gam-1))*(((gam+1)/2)/(gam*x**2 - (gam-1)/2))**(1/(gam-1))    \
          - PtRatio
        Mshock = scipy.optimize.fsolve(func=rootFunc3,x0=2,xtol=1e-12)
        #MpostShock = np.sqrt((1 + ((gam-1)/2)*Mshock**2)/(gam*Mshock**2 -    \
        #  (gam-1)/2))
    elif( pressureRatio < PtRatioSupersonic - deltaPtRatio ):
        status = "overexpanded"
    elif( pressureRatio < PtRatioSupersonic + deltaPtRatio ):
        status = "fully expanded"
    else:
        status = "underexpanded"
        
    return (status,Mshock)
    
#==============================================================================
# Solve for location where M = 1 (location of apparent throat)
#==============================================================================
def findApparentThroat(nozzle,tol,(xInterp,Cf,Tstag,dTstagdx)):
    gam = nozzle.fluid.gam
    relTol = tol["solverApparentThroatLocation"]

    dMdxCoeffFunc = lambda x: -nozzle.wall.geometry.areaGradient(x)/         \
      nozzle.wall.geometry.area(x) + 2*gam*np.interp(x,xInterp,Cf)/          \
      nozzle.wall.geometry.diameter(x) + (1+gam)*np.interp(x,xInterp,        \
      dTstagdx)/(2*np.interp(x,xInterp,Tstag))
      
    # Find sign changes in dMdxCoeffFunc
    xFind = np.linspace(0.,nozzle.wall.geometry.length-1e-4,200.)
    coeffFind = dMdxCoeffFunc(xFind)
    coeffFindSign = np.sign(coeffFind)
    coeffFindSignChange = ((np.roll(coeffFindSign,1) - coeffFindSign)        \
      != 0).astype(int)
    coeffFindSignChange[0] = 0
    signChangeLocations = np.nonzero(coeffFindSignChange)[0]
    
    if( signChangeLocations.size == 0 ):
        throatGuess = nozzle.wall.geometry.findThroat()[0]
    else:
        minInd = np.argmin(nozzle.wall.geometry.area(xFind[signChangeLocations]))
        throatGuess = xFind[signChangeLocations[minInd]]
        
        # Check to make sure each following possible throat is far enough away
        ind = signChangeLocations[minInd]
        for ii in range(minInd+1,signChangeLocations.size):
            currentInd = signChangeLocations[ii]
            dx = xFind[currentInd] - xFind[ind]
            
            #dAdxbar = nozzle.wall.geometry.area(xFind[ind])
            dAdxbarEst = (nozzle.wall.geometry.area(xFind[currentInd]) -     \
              nozzle.wall.geometry.area(xFind[ind]))/dx
            
            dTstagdxbar = np.interp(xFind[currentInd],xInterp,dTstagdx)
            Abar = nozzle.wall.geometry.area(xFind[ind])
            Tstagbar = np.interp(xFind[ind],xInterp,Tstag)
            Cfbar = np.interp(xFind[ind],xInterp,Cf)
            Dbar = nozzle.wall.geometry.diameter(xFind[ind])
            
            RHS = dTstagdxbar*Abar*(1+gam)/(2*Tstagbar) + 2*gam*Abar*Cfbar/Dbar
            
            if( dAdxbarEst <= RHS ):
                throatGuess = xFind[currentInd]
            
        # END OF for ii in range(minInd,signChangeLocations.size)
            
    xApparentThroat = scipy.optimize.fsolve(func=dMdxCoeffFunc,              \
      x0=throatGuess,xtol=relTol)
    
    # Perform some error checking for the apparent throat location
    if( np.isnan(xApparentThroat) or                                         \
      xApparentThroat > nozzle.wall.geometry.length):
        xApparentThroat = nozzle.wall.geometry.length
        print "Mach = 1 at exit\n"
    
    if( nozzle.wall.geometry.diameter(xApparentThroat) >                      \
      nozzle.wall.geometry.diameter(nozzle.wall.geometry.length)):
        xApparentThroat = nozzle.wall.geometry.length

    return xApparentThroat

#==============================================================================
# Define non-ideal quasi-1D equations of motion for post-throat region    
#==============================================================================
def dM2dxPost(x,M2,gam,geo,(xInterp,CfInterp,TstagInterp,dTstagdxInterp)):
    dAdx = geo.areaGradient(x)
    A = geo.area(x)
    D = geo.diameter(x)
    Cf = np.interp(x,xInterp,CfInterp)
    Tstag = np.interp(x,xInterp,TstagInterp)
    dTstagdx = np.interp(x,xInterp,dTstagdxInterp)
    dM2dx = (2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx/A + 2*gam*M2*Cf/D +        \
      (1+gam*M2)*dTstagdx/(2*Tstag))
    return dM2dx
    
def dM2dxPost2(M2,x,gam,geo,(xInterp,CfInterp,TstagInterp,dTstagdxInterp)):
    dAdx = geo.areaGradient(x)
    A = geo.area(x)
    D = geo.diameter(x)
    Cf = np.interp(x,xInterp,CfInterp)
    Tstag = np.interp(x,xInterp,TstagInterp)
    dTstagdx = np.interp(x,xInterp,dTstagdxInterp)
    dM2dx = (2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx/A + 2*gam*M2*Cf/D +        \
      (1+gam*M2)*dTstagdx/(2*Tstag))
    return dM2dx

#==============================================================================
# Integrate subsonic flow through an axial nozzle geometry
#==============================================================================
def integrateSubsonic(nozzle,tol,params):
    if( hasattr(nozzle.inlet, "mach") ): # then integrate using given Mach
        print "Using prescribed inlet Mach number\n"
        M0 = nozzle.inlet.mach
        
#        start = time.time()
#        r = scipy.integrate.ode(dM2dxPost,jac=None).set_integrator('dopri5', \
#          atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
#        r.set_initial_value(M0**2,0.)
#        r.set_f_params(nozzle.fluid.gam,nozzle.wall.geometry,params)
#        dt = nozzle.wall.geometry.length/200.
#        while r.successful() and r.t < nozzle.wall.geometry.length:
#            r.integrate(r.t+dt)
#            #print("%g %g" % (r.t, r.y))
#        M2 = r.integrate(nozzle.wall.geometry.length)
#        end= time.time()
#        print(end-start)
        
#        start = time.time()
        xIntegrate = np.linspace(0.,nozzle.wall.geometry.length,200)
        M2 = scipy.integrate.odeint(dM2dxPost2,M0**2,xIntegrate,             \
          (nozzle.fluid.gam,nozzle.wall.geometry,params),                    \
          atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
#        end = time.time()
#        print(end-start)
        
    else:
        print "else"
    
    return (xIntegrate, M2)

#==============================================================================
# Integrate subsonic, shock post-throat, supersonic flow through an axial 
# nozzle geometry
#==============================================================================
def integrateShock(nozzle,tol,params):
    pass

#==============================================================================
# Integrate supersonic flow through an axial nozzle geometry
#==============================================================================
def integrateSupersonic(nozzle,tol,params):
    pass

#==============================================================================
# Perform quasi-1D area-averaged Navier-Stokes analysis of axisymmetric nozzle.
#% Solve for flow along length of non-ideal nozzle given geometry, inlet
#% stagnation temperature and pressure, and freestream temperature and
#% pressure. Iterate for Cf and stagnation temperature. An ODE for M^2 is 
#% solved given A, Cf, and Tstag. Pstag is found from mass conservation. T 
#% and P are found from def'n of stag. temp. Density rho is found from ideal 
#% gas law. 
#%
#% Returns M, density, pressure P, temperature T, stagnation 
#% temp. Tstag, stagnation pressure Pstag, velocity U, Re, internal heat
#% transfer coefficient hf, friction coefficient Cf, interior wall temp. Tw,
#% exterior wall temp. Text, and approximate stress along length of nozzle.
#==============================================================================
def analysis(nozzle,tol):
    
    # Initialize
    xApparentThroat = nozzle.wall.geometry.findMinimumRadius()[0]
    
    # Determine state of nozzle assuming ideal conditions
    pressureRatio = nozzle.inlet.Pstag/nozzle.environment.P
    TstagThroat = nozzle.inlet.Tstag
    PstagThroat = nozzle.inlet.Pstag
    TstagExit = nozzle.inlet.Tstag
    PstagExit = nozzle.inlet.Pstag    
    (status,shock) = nozzleState(nozzle,pressureRatio,PstagThroat,           \
      TstagThroat,PstagExit,TstagExit)
    
    # Initialize loop variables
    xInterp = np.array(([0., nozzle.wall.geometry.length]))
    Cf = np.array(([0.004, 0.004]))
    Tstag = np.array(([nozzle.inlet.Tstag, nozzle.inlet.Tstag]))
    dTstagdx = np.array(([-6., -6.]))
    xPositionOld = np.array(([0., nozzle.wall.geometry.length]))
    
    maxIterations = 10 # max number of iterations to solve for Cf and Tstag
    counter = 0 # used to count b/w number of iterations
    tolerance = tol["exitTempPercentError"] # tolerance for % error in
                # exit static temperature between iterations
    Texit_old = 0 # save previous static temperature
    
    while( 1 ):
        
        counter += 1
        
        # ODE solver options
        odeSolverRelTol = tol["solverRelTol"]
        odeSolverAbsTol = tol["solverAbsTol"]
        #odeEvents = define eventsFcn
        
        params = (xInterp,Cf,Tstag,dTstagdx)
        
        # Find where M = 1
        xApparentThroat = findApparentThroat(nozzle,tol,params)
        
        if( status == "no flow" ):
            raise UserWarning("Prescribed inputs result in flow reversal     \
              in nozzle")
        elif( status == "subsonic" ):
            (xPosition,M2) = integrateSubsonic(nozzle,tol,params)
        elif( status == "shock" ):
            (xPosition,M2) == integrateShock(nozzle,tol,params)
        else: # supersonic flow
            (xPosition,M2) == integrateSupersonic(nozzle,tol,params)
            
        print xPosition
        print M2
        
        
        
        break
    
    # END OF while( ~converged )
    
# END OF analysis(nozzle,tol)
    