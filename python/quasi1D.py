# -*- coding: utf-8 -*-
"""
Perform quasi-1D area-averaged Navier-Stokes analysis on axisymmetric nozzle

Rick Fenrich 6/28/16
"""

import numpy as np
import scipy.optimize
import scipy.integrate

import lifetime
import geometry

#==============================================================================
# Sutherland's Law of dynamic viscosity of air
#==============================================================================
def dynamicViscosity(T):
    mu = 1.716e-5*(T/273.15)**1.5*(273.15 + 110.4)/(T + 110.4)
    return mu
   
#==============================================================================
# Area-Mach function from 1-D mass conservation equations
#==============================================================================
def areaMachFunc(g,M):
    a = ((g+1)/2)**((g+1)/(2*(g-1)))*M/(1+(g-1)*M**2/2)**((g+1)/(2*(g-1)))
    return a

#==============================================================================
# Mass flow-rate
#==============================================================================
def massFlowRate(fluid,Pstag,Area,Tstag,M):
    gam = fluid.gam
    R = fluid.R;
    mdot = (gam/((gam+1)/2)**((gam+1)/(2*(gam-1))))*Pstag*Area*              \
      areaMachFunc(gam,M)/np.sqrt(gam*R*Tstag)
    return mdot

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
        throatGuess = nozzle.wall.geometry.findMinimumRadius()[0]
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
      x0=throatGuess,xtol=relTol)[0]
    
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
# Define non-ideal quasi-1D equations of motion for forward integration    
#==============================================================================
def dM2dxForward(x,M2,gam,geo,(xInterp,CfInterp,TstagInterp,dTstagdxInterp)):
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
# Define non-ideal quasi-1D equations of motion for backward integration    
#==============================================================================    
def dM2dxBackward(x,M2,xRef,gam,geo,(xInterp,CfInterp,TstagInterp,dTstagdxInterp)):
    dAdx = geo.areaGradient(xRef-x)
    A = geo.area(xRef-x)
    D = geo.diameter(xRef-x)
    Cf = np.interp(xRef-x,xInterp,CfInterp)
    Tstag = np.interp(xRef-x,xInterp,TstagInterp)
    dTstagdx = np.interp(xRef-x,xInterp,dTstagdxInterp)
    dM2dx = -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx/A + 2*gam*M2*Cf/D +       \
      (1+gam*M2)*dTstagdx/(2*Tstag))
    return dM2dx
    
#==============================================================================
# Check ODE integration for problems such as hitting a singularity
#==============================================================================
#def checkODEintegration(x,y):
#    r = 2. # maximum allowable multiplicative factor for slopes between 2
#            # adjacent steps
#    tol = 0.01 # tolerance to determine whether y is close enough to zero
#    
#    # calculate slopes for each step
#    m = (y - np.roll(y,1))/(x - np.roll(x,1))
#    m[0] = m[1] # assume problem doesn't occur at first step
#
#    v1 = [ii for ii in range(1,m.size) if abs(m[ii]) > r*abs(m[ii-1]) and (m[ii-1] > tol or m[ii] > tol)]
#    
#    if( len(v1) != 0 ):
#        minIndex = v1[0] - 1 # last good value before solver went crazy
#    else:
#        minIndex = -1
#
#    return minIndex
    
#==============================================================================
# Integrate ODE built with scipy.integrate.ode and perform simple event check.
# Provide solution at equally spaced points. Integration will end when event is
# detected.
#==============================================================================
def integrateODEwithEvents(ode,dt,tfinal,ycrit,direction):
    
    if( isinstance(dt,float) ): # uniform spacing
        if( direction == "b" ):
            x = np.linspace(tfinal,ode.t,np.round((ode.t-tfinal)/dt))
        else:
            x = np.linspace(ode.t,tfinal,np.round((tfinal-ode.t)/dt))
    else:
        raise TypeError(("integration for y(t) at uniform intervals not ",
                        "enabled, dt must be a float"))
        # if non-uniform x-spacing is desired, x must be flipped for backwards
        # integration
        # for backwards integration: x = abs(x - max(x)) + min(x)
    
    ii = 0
    y = np.zeros(x.size)
    delta = 1e-12 # gives distance away from event when solver terminates
    
    if( direction == "b" ): # backward integration
    
        # translate into a forward integration problem
        tTemp = ode.t
        ode.set_initial_value(ode.y,tfinal)
        tfinal = tTemp
    
        if( ode.y > ycrit ):
            
            while ode.successful() and ode.t < tfinal - dt/2:
                
                yTemp = ode.integrate(ode.t+dt)
                
                if( yTemp < ycrit + delta ):
                    y[x.size-ii-1] = yTemp
                    break # critical y reached
                    
                y[x.size-ii-1] = yTemp
                ii +=1
                    
        elif( ode.y <= ycrit ):
        
            while ode.successful() and ode.t < tfinal - dt/2:
                
                yTemp = ode.integrate(ode.t+dt)
                
                if( yTemp > ycrit - delta ):
                    y[x.size-ii-1] = yTemp
                    break # critical y reached
                
                y[x.size-ii-1] = yTemp
                ii +=1       
    
    else: # forward integration
    
        if( ode.y > ycrit ):
            
            while ode.successful() and ode.t < tfinal - dt/2:
                
                yTemp = ode.integrate(ode.t+dt)
                
                if( yTemp < ycrit + delta ):
                    y[ii] = yTemp
                    break # critical y reached
                    
                y[ii] = yTemp
                ii +=1
                    
        elif( ode.y <= ycrit ):
        
            while ode.successful() and ode.t < tfinal - dt/2:
                
                yTemp = ode.integrate(ode.t+dt)
                
                if( yTemp > ycrit - delta ):
                    y[ii] = yTemp
                    break # critical y reached
                
                y[ii] = yTemp
                ii +=1  
    
    return (x,y,ii)
    
#==============================================================================
# Integrate subsonic flow through an axial nozzle geometry
#==============================================================================
def integrateSubsonic(nozzle,tol,params,xThroat,nPartitions):
    
    # Use inlet Mach number provided by user, if available
    if( hasattr(nozzle.inlet, "mach") ):
        print "Using prescribed inlet Mach number\n"
        M0 = nozzle.inlet.mach
        
        f = scipy.integrate.ode(dM2dxForward,jac=None).set_integrator(       \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        f.set_initial_value(M0**2,0.)
        f.set_f_params(nozzle.fluid.gam,nozzle.wall.geometry,params)
        dt = nozzle.wall.geometry.length/nPartitions
        tfinal = nozzle.wall.geometry.length
        (xIntegrate,M2,eventIndex) = integrateODEwithEvents(f,dt,tfinal,1.,"f")
        
        if( eventIndex != xIntegrate.size ):
            raise ValueError(("Integration terminated early: prescribed ",
                              "inlet Mach number is too large"))
    
    # Else if nozzle converges only, assume choked flow at the exit    
    elif( nozzle.wall.geometry.length - xThroat < 1e-12 ):
        M0 = 0.9999 # start integration from this Mach number
        
        b = scipy.integrate.ode(dM2dxBackward,jac=None).set_integrator(      \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        b.set_initial_value(M0**2,nozzle.wall.geometry.length)
        b.set_f_params(nozzle.wall.geometry.length,nozzle.fluid.gam,         \
          nozzle.wall.geometry,params)
        dt = nozzle.wall.geometry.length/nPartitions
        tfinal = 0.
        (xIntegrate,M2,eventIndex) = integrateODEwithEvents(b,dt,tfinal,1.,"b")
                
        if( eventIndex != xIntegrate.size ):
            raise RuntimeError(("Integration terminated early while ",
                               "integrating backwards from the exit"))
    
    # Else, assume nozzle is choked at throat, integrate forwards & backwards
    else:
        UpperM = 0.9999 # start integration at this Mach number for aft portion
        LowerM = 0.9999 # start integ. at this Mach number for fore portion
        dx = 1e-5 # 1e-5 for 0.9999 to 1.0001 or 1e-4 for 0.999 to 1.001
        
        # Integrate forward from throat
        f = scipy.integrate.ode(dM2dxForward,jac=None).set_integrator(       \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        f.set_initial_value(UpperM**2,xThroat+dx/2)
        f.set_f_params(nozzle.fluid.gam,nozzle.wall.geometry,params)
        nP = np.round((1-xThroat/nozzle.wall.geometry.length)*nPartitions)
        dt = (nozzle.wall.geometry.length-xThroat-dx/2)/nP
        tfinal = nozzle.wall.geometry.length
        (xF,M2F,eventIndex) = integrateODEwithEvents(f,dt,tfinal,1.,"f")
        
        if( eventIndex != xF.size ):
            raise RuntimeError(("Integration terminated early while ",
                               "integrating forwards from the throat"))
        
        # Integrate backward from throat
        b = scipy.integrate.ode(dM2dxBackward,jac=None).set_integrator(      \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        b.set_initial_value(LowerM**2,xThroat-dx/2)
        b.set_f_params(xThroat,nozzle.fluid.gam,nozzle.wall.geometry,params)
        nP = np.round(xThroat/nozzle.wall.geometry.length*nPartitions)
        dt = (xThroat-dx/2)/nP
        tfinal = 0.
        (xB,M2B,eventIndex) = integrateODEwithEvents(b,dt,tfinal,1.,"b")
                
        if( eventIndex != xB.size ):
            raise RuntimeError(("Integration terminated early while ",
                               "integrating backwards from throat"))
    
        xIntegrate = np.concatenate((xB,xF))
        M2 = np.concatenate((M2B,M2F))
    
    return (xIntegrate, M2)

#==============================================================================
# Integrate subsonic, shock post-throat, supersonic flow through an axial 
# nozzle geometry
#==============================================================================
def integrateShock(nozzle,tol,params,xThroat,nPartitions):
    raise RuntimeError("Shock in nozzle case not implemented yet")
    return -1

#==============================================================================
# Integrate supersonic flow through an axial nozzle geometry
#==============================================================================
def integrateSupersonic(nozzle,tol,params,xThroat,nPartitions):
        
    # If nozzle converges only, assume choked flow at the exit    
    if( nozzle.wall.geometry.length - xThroat < 1e-12 ):
        M0 = 0.9999 # start integration from this Mach number
        
        b = scipy.integrate.ode(dM2dxBackward,jac=None).set_integrator(      \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        b.set_initial_value(M0**2,nozzle.wall.geometry.length)
        b.set_f_params(xThroat,nozzle.wall.geometry.length,nozzle.fluid.gam,         \
          nozzle.wall.geometry,params)
        dt = nozzle.wall.geometry.length/nPartitions
        tfinal = 0.
        (xIntegrate,M2,eventIndex) = integrateODEwithEvents(b,dt,tfinal,1.,"b")
                
        if( eventIndex != xIntegrate.size ):
            raise RuntimeError(("Integration terminated early while ",
                               "integrating backwards from the exit"))

    # Else, assume nozzle is choked at throat, integrate forwards & backwards
    else:
        UpperM = 1.0001 # start integration at this Mach number for aft portion
        LowerM = 0.9999 # start integ. at this Mach number for fore portion
        dx = 1e-4 # 1e-5 for 0.9999 to 1.0001 or 1e-4 for 0.999 to 1.001
        
        # Integrate forward from throat
        f = scipy.integrate.ode(dM2dxForward,jac=None).set_integrator(       \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        f.set_initial_value(UpperM**2,xThroat+dx/2)
        f.set_f_params(nozzle.fluid.gam,nozzle.wall.geometry,params)
        nP = np.round((1-xThroat/nozzle.wall.geometry.length)*nPartitions)
        dt = (nozzle.wall.geometry.length-xThroat-dx/2)/nP
        tfinal = nozzle.wall.geometry.length
        (xF,M2F,eventIndex) = integrateODEwithEvents(f,dt,tfinal,1.,"f")
        
        if( eventIndex != xF.size ):
            raise RuntimeError(("Integration terminated early while ",
                               "integrating forwards from the throat"))
        
        # Integrate backward from throat
        b = scipy.integrate.ode(dM2dxBackward,jac=None).set_integrator(      \
          'dopri5',atol=tol["solverAbsTol"],rtol=tol["solverRelTol"])
        b.set_initial_value(LowerM**2,xThroat-dx/2)
        b.set_f_params(xThroat,nozzle.fluid.gam,nozzle.wall.geometry,params)
        nP = np.round(xThroat/nozzle.wall.geometry.length*nPartitions)
        dt = (xThroat-dx/2)/nP
        tfinal = 0.
        (xB,M2B,eventIndex) = integrateODEwithEvents(b,dt,tfinal,1.,"b")
                
        if( eventIndex != xB.size ):
            raise RuntimeError(("Integration terminated early while ",
                                "integrating backwards from throat"))
    
        xIntegrate = np.concatenate((xB,xF))
        M2 = np.concatenate((M2B,M2F))
    
    return (xIntegrate, M2)

#==============================================================================
# Cumulative trapezoidal integration using linear interpolation to return n 
# results for a vector of size n    
#==============================================================================
def integrateTrapezoidal(y,x):
        integrand = np.empty(x.size+1)
        dx = np.empty(x.size+1)
        
        dx0 = x[1] - x[0]
        dxN = x[-1] - x[-2]
        
        integrand[0] = y[0] - (y[1] - y[0])/dx0*(x[0] - (x[0] - dx0/2))
        integrand[1:-1] = y[1:] - (y[1:] - y[:-1])/(x[1:] - x[:-1])*         \
          (x[1:] - (x[1:] + x[:-1])/2)
        integrand[-1] = y[-1] + (y[-1] - y[-1])/dxN*((x[-1] + dxN/2) - x[-1])
        
        dx[0] = x[0] - dx0/2
        dx[1:-1] = (x[1:] + x[:-1])/2
        dx[-1] = x[-1] + dxN/2
        
        integral = scipy.integrate.cumtrapz(integrand,dx)
        
        return (integral)

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
    gam = nozzle.fluid.gam
    xApparentThroat = nozzle.wall.geometry.findMinimumRadius()[0]
    nOdeIntegrationSteps = 1000;
    
    # Determine state of nozzle assuming ideal conditions
    pressureRatio = nozzle.inlet.Pstag/nozzle.environment.P
    TstagThroat = nozzle.inlet.Tstag
    PstagThroat = nozzle.inlet.Pstag
    TstagExit = nozzle.inlet.Tstag
    PstagExit = nozzle.inlet.Pstag    
    (status,shock) = nozzleState(nozzle,pressureRatio,PstagThroat,           \
      TstagThroat,PstagExit,TstagExit)
    
    # Initialize loop variables
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
        
        # Parameters passed to functions called by ODE        
        params = (xPositionOld,Cf,Tstag,dTstagdx)
        
        # Find where M = 1
        xApparentThroat = findApparentThroat(nozzle,tol,params)
        
        if( status == "no flow" ):
            raise UserWarning(("Prescribed inputs result in flow reversal ",
                               "in nozzle"))
        elif( status == "subsonic" ):
            (xPosition,M2) = integrateSubsonic(nozzle,tol,params,            \
              xApparentThroat,nOdeIntegrationSteps)
        elif( status == "shock" ):
            (xPosition,M2) = integrateShock(nozzle,tol,params,               \
              xApparentThroat,nOdeIntegrationSteps)
        else: # supersonic flow
            (xPosition,M2) = integrateSupersonic(nozzle,tol,params,          \
              xApparentThroat,nOdeIntegrationSteps)
              
        # Check output
        if( np.isnan(M2.any()) or M2.any() < 0. or np.isinf(M2.any()) ):
            raise RuntimeError("Unrealistic Mach number calculated")
            
        # Calculate geometric properties
        D = nozzle.wall.geometry.diameter(xPosition)
        A = nozzle.wall.geometry.area(xPosition)
        #dAdx = nozzle.wall.geometry.dAdx(xPosition)
        t = nozzle.wall.thickness.radius(xPosition)
            
        # Calculate other 1D flow properties
        M = np.sqrt(M2)
        Tstag = np.interp(xPosition,xPositionOld,Tstag)
        T = Tstag/(1. + (gam-1.)*M2/2.) # static temp. from stag. temp. def.
        Pstag = nozzle.inlet.Pstag*(A[0]/A)*(areaMachFunc(gam,M[0])/         \
          areaMachFunc(gam,M))*np.sqrt(Tstag/Tstag[0]) # from mass conserv.
        P = Pstag/(1. + (gam-1.)*M2/2.)**(gam/(gam-1)) # from stag. press. def.
        density = P/(nozzle.fluid.R*T)
        U = M*np.sqrt(gam*nozzle.fluid.R*T) # velocity
        Re = density*U*D/dynamicViscosity(T) # Reynolds number from definition
        Cf = np.interp(xPosition,xPositionOld,Cf) # friction coefficient
            
        # Recalculate friction and heat
        # heat transfer coefficient to interior nozzle wall, estimated using
        # Chilton-Colburn analogy
        hf = nozzle.fluid.Pr(T)**(2./3.)*density*nozzle.fluid.Cp(T)*U*Cf/2
        
        # Redefine stagnation temperature distribution
        TstagXIntegrand = 4./(nozzle.fluid.Cp(T)*density*U*D*(1./hf +        \
          t/nozzle.wall.material.k + 1./nozzle.environment.hInf))   
        TstagXIntegral = integrateTrapezoidal(TstagXIntegrand,xPosition)
        Tstag = nozzle.environment.T*(1. - np.exp(-TstagXIntegral)) +        \
          nozzle.inlet.Tstag*np.exp(-TstagXIntegral)
        dTstagdx = (nozzle.environment.T - Tstag)*4./(nozzle.fluid.Cp(T)*    \
          density*U*D*(1./hf + t/nozzle.wall.material.k +                    \
          1./nozzle.environment.hInf))
          
        # Estimate interior wall temperature
        Qw = nozzle.fluid.Cp(T)*density*U*D*dTstagdx/4.
        Tinside = Tstag + Qw/hf # interior wall temperature
        #recoveryFactor = (Tinside/T - 1)/((gam-1)*M2/2)
        
        # Estimate exterior wall temperature
        Toutside = nozzle.environment.T - Qw/nozzle.environment.hInf
    
        # Redefine friction coefficient distribution (Sommer & Short's method)
        TPrimeRatio = 1. + 0.035*M2 + 0.45*(Tinside/T - 1.)
        RePrimeRatio = 1./(TPrimeRatio*(TPrimeRatio)**1.5*(1. + 110.4/T)/    \
          (TPrimeRatio + 110.4/T))
        CfIncomp = 0.074/Re**0.2
        Cf = CfIncomp/TPrimeRatio/RePrimeRatio**0.2
        
        # Save old solution x position for next iteration
        xPositionOld = xPosition    
        
        # Recalculate nozzle status
        TstagThroat = np.interp(xApparentThroat,xPosition,Tstag)
        PstagThroat = np.interp(xApparentThroat,xPosition,Pstag)
        TstagExit = Tstag[-1]
        if( status == "shock" ):
            # approximate calculation
            shockIndex = 10000 
            pressureRatio = Pstag[shockIndex-1]/nozzle.environment.P
            PstagExit = Pstag[shockIndex-1] # temp. fix
            (status,shock) = nozzleState(nozzle,pressureRatio,PstagThroat,   \
              TstagThroat,PstagExit,TstagExit)
            PstagExit = Pstag[-1]
        else: # no shock in nozzle; assume no Pstag loss
            pressureRatio = Pstag[-1]/nozzle.environment.P
            PstagExit = Pstag[-1]
            (status,shock) = nozzleState(nozzle,pressureRatio,PstagThroat,   \
              TstagThroat,PstagExit,TstagExit)
        
        #print "%i\n" % counter
        
        if( counter >= maxIterations ):
            print "Iteration limit for quasi-1D heat xfer & friction reached\n"
            break
        
        # Check tolerance on static temperature at nozzle exit
        percentError = abs(T[-1] - Texit_old)/T[-1]
        #print "Percent error: %e\n" % (percentError*100)
        Texit_old = T[-1]
        if( percentError < tolerance ):
            #print "%i iterations to converge quasi-1D heat xfer & friction \
#calcs\n" % counter    
            break
    
    # END OF while( ~converged )
    
    dAdx = nozzle.wall.geometry.areaGradient(xPosition)
    maxSlope = max(dAdx/np.pi/D)
    minSlope = min(dAdx/np.pi/D)
    
    # Calculate mass flow rate
    mdot = massFlowRate(nozzle.fluid,Pstag,A,Tstag,M)
    
    # Calculate thrust
    exitAngle = np.arctan2(dAdx[-1],D[-1])
    divergenceFactor = (1. + np.cos(exitAngle))/2.
    netThrust = divergenceFactor*mdot[0]*(U[-1] - nozzle.mission.mach*       \
      nozzle.environment.c) + (P[-1] - nozzle.environment.P)*A[-1]
    #grossThrust = divergenceFactor*mdot[0]*U + (P[-1] -                      \
    #  nozzle.environment.P)*A[-1]
      
    # Calculate stresses
    # Stresses calculated assuming cylinder; nozzle length not constrained in 
    # thermal expansion
    stressHoop = P*D/(2.*t)
    # Thermal stresses calculated assuming steady-state, give max tens. stress
    ri = D/2. # inner radius
    ro = D/2. + t # outer radius
    stressThermalRadial = nozzle.wall.material.E*                            \
      nozzle.wall.material.alpha*(Tinside-Toutside)/                         \
      (2.*(1.-nozzle.wall.material.v))*(1./np.log(ro/ri))*(1. - 2.*ri**2./   \
      (ro**2. - ri**2.)*np.log(ro/ri))
    stressThermalTangential = stressThermalRadial
    # Estimate vonMises, even though not really valid for composites
    #stressVonMises = np.sqrt( (stressHoop+stressThermalTangential)**2 -      \
    #   stressHoop*stressThermalRadial + stressThermalRadial**2 )
    stressMaxPrincipal = stressHoop + stressThermalTangential
    stressPrincipal = (stressMaxPrincipal, stressThermalRadial,              \
      np.zeros(xPosition.size))
      
    # Calculate cycles to failure, Nf
    Nf = lifetime.estimate(Tinside,stressMaxPrincipal,1)
    
    # Calculate volume of nozzle material (approximately using trap. integ.)
    volume = geometry.wallVolume(nozzle.wall.geometry,nozzle.wall.thickness)
    
    # Assign all data for output
    #flowDimension = 1 # 1-dimensional flow field
    
    flowTuple = (M, U, density, P, Pstag, T, Tstag, Re)
    heatTuple = (Tinside, Toutside, Cf, hf, Qw)
    geoTuple = (D, A, t, dAdx, minSlope, maxSlope)
    performanceTuple= (volume, netThrust, Nf, mdot, Pstag[-1]/Pstag[0],      \
      Tstag[-1]/Tstag[0], status)
    stressTuple = (stressHoop, stressThermalRadial, stressThermalTangential, \
      stressMaxPrincipal, stressPrincipal)
      
    return (xPosition, flowTuple, heatTuple, geoTuple, stressTuple,          \
      performanceTuple)
    
# END OF analysis(nozzle,tol)
    