# -*- coding: utf-8 -*-
"""
Geometry module for axisymmetric nozzle. Currently only a B-spline geometry is 
implemented.

Rick Fenrich 6/28/16
"""

import numpy as np 
import scipy.optimize
import scipy.integrate   
import geometryC

class Bspline():
    def __init__(self, coefs): # assumes 3rd degree B-spline
        self.type = "B-spline"
        self.coefs = coefs
        self.knots = np.hstack(([np.zeros(4), np.arange(1.,coefs.size/2-3),  \
          np.ones(4)*(coefs.size/2-3)])) # calculate here
        self.degree = self.knots.size - self.coefs.size/2 - 1
        self.length = coefs[0,-1]
        self.inletRadius = coefs[1,0]
        
    def findMinimumRadius(self):
        xSeg = np.zeros(self.knots.size)
        ySeg = np.zeros(self.knots.size)
        for ii in range(0,self.knots.size):
            (xTemp,yTemp,temp1,temp2) = uMap3(self.knots[ii],self)
            xSeg[ii] = xTemp
            ySeg[ii] = yTemp
        yMinKnotIndex = np.argmin(ySeg)
        minFunc = lambda x: self.diameter(x).item()
        DMinOld = 1e12
        for ii in range(max(yMinKnotIndex-2,0),min(yMinKnotIndex+3,          \
          self.knots.size-1)):
            lowerBound = xSeg[ii].item()
            upperBound = xSeg[ii+1].item()
            xMin = scipy.optimize.fminbound(minFunc,lowerBound,upperBound)  
            DMin = minFunc(xMin)
            if( DMin < DMinOld ):
                DMinOld = DMin
                xMinOld = xMin
        self.xThroat = xMinOld
        self.yThroat = DMinOld/2
        self.Ainlet2Athroat = (self.inletRadius)**2/self.yThroat**2
        self.Aexit2Athroat = (self.coefs[1,-1])**2/self.yThroat**2
        return (self.xThroat, self.yThroat)
        
    def radius(self, x): # r
        #y = bSplineGeometry(x,self)[0] # Python version (slower)
        y = bSplineGeometryC(x,self)[0] # uses dynamic C library
        return y     
        
    def diameter(self, x): # D
        #y = bSplineGeometry(x,self)[0] # Python version (slower)
        y = bSplineGeometryC(x,self)[0] # uses dynamic C library
        return y*2
        
    def area(self, x): # A
        #y = bSplineGeometry(x,self)[0] # Python version (slower)
        y = bSplineGeometryC(x,self)[0] # uses dynamic C library
        return np.pi*y**2
        
    def areaGradient(self, x): # dAdx
        #(y, dydx) = bSplineGeometry(x,self) # Python version (slower)
        (y, dydx) = bSplineGeometryC(x,self) # uses dynamic C library
        return 2*np.pi*y*dydx
        
class PiecewiseLinear:
    def __init__(self,nodes):
        self.type = "piecewise-linear"
        self.nodes = nodes
        self.length = nodes[0,-1]
        self.inletRadius = nodes[1,0]
        
    def findMinimumRadius(self):
        ii = np.argmin(self.nodes[1,:])
        self.xThroat = self.nodes[0,ii]
        self.yThroat = self.nodes[1,ii]
        self.Ainlet2Athroat = (self.inletRadius)**2/self.yThroat**2
        self.Aexit2Athroat = (self.nodes[1,-1])**2/self.yThroat**2
        return (self.xThroat, self.yThroat)
        
    def radius(self, x): # r
        y = np.interp(x,self.nodes[0,:],self.nodes[1,:])
        return y
        
    def diameter(self, x): # D
        y = np.interp(x,self.nodes[0,:],self.nodes[1,:])
        return 2*y
        
    def area(self, x): # A
        y = np.interp(x,self.nodes[0,:],self.nodes[1,:])
        return np.pi*y**2
        
    def areaGradient(self, x): # dAdx
        y = np.interp(x,self.nodes[0,:],self.nodes[1,:])
        if( isinstance(x,float) ):
            upperIndex = find(x,self.nodes[0,:])
            if( upperIndex == self.nodes.size/2 ):
                upperIndex = upperIndex - 1
            lowerIndex = upperIndex - 1
            dydx = (self.nodes[1,upperIndex] - self.nodes[1,lowerIndex])/    \
              (self.nodes[0,upperIndex] - self.nodes[0,lowerIndex])
        else: # x is an array
            dydx = np.zeros(x.size)
            for ii in range(0,x.size):
                upperIndex = find(x[ii],self.nodes[0,:])
                if( upperIndex == self.nodes.size/2 ):
                    upperIndex = upperIndex - 1
                lowerIndex = upperIndex - 1
                dydx[ii] = (self.nodes[1,upperIndex] - 
                  self.nodes[1,lowerIndex])/(self.nodes[0,upperIndex]        \
                  - self.nodes[0,lowerIndex])
            
        return 2*np.pi*y*dydx

#==============================================================================
# Find first 1-based index where scalar xFind < xVec[ii]
#==============================================================================
def find(xFind,xVec):
    counter = 1
    for ii in range(1,xVec.size):
        if(xFind < xVec[ii]):
            break # assume xVec is in ascending order
        counter += 1
    return counter

#==============================================================================
# Calculate x, y, dxdu, and dydu given u for 3rd degree B-spline
# Currently implemented only for scalar inputs of u
#==============================================================================
def uMap3(u,bSpline):
    x = 0.
    y = 0.
    dxdu = 0.
    dydu = 0.
    
    # Calculate the highest knot index that gives a value of u below the given
    # value (1-based index)
    hh = find(u,bSpline.knots)
    if( hh == bSpline.knots.size ): # at end of knots vector
        hh = bSpline.knots.size - 4
        
    nn = 0
    if( hh == 1 ):
        nn = 1
    elif( hh == 2 ):
        nn = 2
    elif( hh == 3 ):
        nn = 3
    else:
        nn = 4
        
    ii = 0
    while( ii > -nn ): # i.e. for each contributing basis
        jj = hh + ii - 1 # subtracted 2 for C++ compatibility
        
        # Redefine k1 through k5 here
        k1 = bSpline.knots[jj]
        k2 = bSpline.knots[jj+1]
        k3 = bSpline.knots[jj+2]
        k4 = bSpline.knots[jj+3]
        k5 = bSpline.knots[jj+4]
        
        if( ii == 0 ): # calculate basis N1
            if( abs(k1-k2) <= 1e-12 ):
                Ncurrent = 0.
                dNducurrent = 0.
            else:
                Ncurrent = (u-k1)/(k4-k1)*(u-k1)/(k3-k1)*(u-k1)/(k2-k1)
                dNducurrent = -(3*pow(k1 - u,2))/((k1 - k2)*(k1 - k3)*       \
                  (k1 - k4))
        elif( ii == -1 ): # calculate basis N2
            if( abs(k2-k3) <= 1e-12):
                Ncurrent = 0.
                dNducurrent = 0.
            else:
                Ncurrent = (u-k1)/(k4-k1)*((u-k1)/(k3-k1)*(k3-u)/(k3-k2) +   \
                  (k4-u)/(k4-k2)*(u-k2)/(k3-k2)) + (k5-u)/(k5-k2)*(u-k2)     \
                  /(k4-k2)*(u-k2)/(k3-k2)
                dNducurrent = (((k1 - u)*(k3 - u))/((k1 - k3)*(k2 - k3)) +   \
                  ((k2 - u)*(k4 - u))/((k2 - k3)*(k2 - k4)))/(k1 - k4) +     \
                  pow(k2 - u,2)/((k2 - k3)*(k2 - k4)*(k2 - k5)) + ((k5 - u)  \
                  *(2*k2 - 2*u))/((k2 - k3)*(k2 - k4)*(k2 - k5)) +           \
                  (2*(k1 - u)*(k1*k2 - k3*k4 - k1*u - k2*u + k3*u + k4*u))/  \
                  ((k1 - k3)*(k1 - k4)*(k2 - k3)*(k2 - k4))
        elif( ii == -2 ): # calculate basis N3
            if( abs(k3-k4) <= 1e-12 ):
                Ncurrent = 0.
                dNducurrent = 0.
            else:
                Ncurrent = (u-k1)/(k4-k1)*(k4-u)/(k4-k2)*(k4-u)/(k4-k3) +    \
                  (k5-u)/(k5-k2)*((u-k2)/(k4-k2)*(k4-u)/(k4-k3) + (k5-u)/    \
                  (k5-k3)*(u-k3)/(k4-k3))
                dNducurrent = - (((k2 - u)*(k4 - u))/((k2 - k4)*(k3 - k4)) + \
                  ((k3 - u)*(k5 - u))/((k3 - k4)*(k3 - k5)))/(k2 - k5) -     \
                  pow(k4 - u,2)/((k1 - k4)*(k2 - k4)*(k3 - k4)) -            \
                  ((k1 - u)*(2*k4 - 2*u))/((k1 - k4)*(k2 - k4)*(k3 - k4)) -  \
                  (2*(k5 - u)*(k2*k3 - k4*k5 - k2*u - k3*u + k4*u + k5*u))/  \
                  ((k2 - k4)*(k2 - k5)*(k3 - k4)*(k3 - k5))
        else: # calculate basis N4
            if( abs(k4-k5) <= 1e-12 ):
                Ncurrent = 0.
                dNducurrent = 0.
            else:            
                Ncurrent = (k5-u)/(k5-k2)*(k5-u)/(k5-k3)*(k5-u)/(k5-k4)
                dNducurrent = (3*pow(k5 - u,2))/((k2 - k5)*(k3 - k5)*        \
                  (k4 - k5))
        
        x += bSpline.coefs[0,jj]*Ncurrent
        y += bSpline.coefs[1,jj]*Ncurrent
        dxdu += bSpline.coefs[0,jj]*dNducurrent
        dydu += bSpline.coefs[1,jj]*dNducurrent
        
        ii -= 1
        
        if( ii < -4):
            break
        
    # END OF while( ii > -nn )
        
    return (x,y,dxdu,dydu)
    
# END OF uMap3

#==============================================================================
# Return y given x for a 3rd degree B-spline
#==============================================================================
def bSplineGeometry(x,bSpline):
    
    if( isinstance(x,float) ):
        if( x > (bSpline.coefs[0][-1]) ):
            x = bSpline.coefs[0][-1]
        elif( x < (bSpline.coefs[0][0]) ):
            x = bSpline.coefs[0][0]
            #raise ValueError("x is outside bounds of B-spline")
    
    if( isinstance(x,np.ndarray) ):
        nx = x.size # number of x
    else: # convert to numpy array for calculations
        x = np.array(([x]))
        nx = 1
        
    k = bSpline.knots.size # number of knots
    
    y = np.zeros(nx)
    dydx = np.zeros(nx)
    
    # Determine x-value at breaks
    xKnot = np.empty(k)
    for ii in range(0,k):
        xKnot[ii] = uMap3(bSpline.knots[ii],bSpline)[0]
    xKnot[k-4] += 1e-6 # so finding upper bound on u works (assumes same 4
                       # knots at end of knot vector)
    
    tolerance = 1e-6 # tolerance for Newton solver
    
    for ii in range(0,nx):
        
        # Determine lower and upper bounds on u
        seg = find(x[ii],xKnot)
        uLower = bSpline.knots[seg-1]
        uUpper = bSpline.knots[seg]
        
        # Pick a guess for u (a linear interpolation)
        u = (uLower + uUpper)/2
        
        # Calculate x and dxdu corresponding to u
        (xEst, temp1, dxduEst, temp2) = uMap3(u,bSpline)
        
        # Perform 1 Newton iteration
        if( dxduEst < tolerance ):
            uNew = 0.
        else:
            uNew = u - (xEst - x[ii])/dxduEst
            
        # Perform remaining Newton iterations
        counter = 0
        errorMeasure = abs((uNew - u)/uNew)
        while( errorMeasure > tolerance ):
            u = uNew
            
            (xEst, temp1, dxduEst, temp2) = uMap3(u,bSpline)
            
            if( dxduEst < 1e-12 ):
                uNew = u
            else:
                uNew = u - (xEst - x[ii])/dxduEst
                
            counter += 1
            
            if( counter > 20 ):
                break
            
            errorMeasure = abs((uNew-u)/uNew)
            
        # END OF while( errorMeasure > tolerance )
            
        u = uNew
        
        (xTemp, yTemp, dxduTemp, dyduTemp) = uMap3(u,bSpline)
        
        y[ii] = yTemp
        
        if( dxduTemp < tolerance ):
            dydx[ii] = 0.
        else:
            dydx[ii] = dyduTemp/dxduTemp
            
    # END OF for ii in range(0,nx)
            
    return (y, dydx)

# END OF bSplineGeometry
        
#==============================================================================
# Return y given x for a 3rd degree B-spline
#==============================================================================
def bSplineGeometryC(x,bSpline):
    
    if( isinstance(x,float) ):
        if( x > (bSpline.coefs[0][-1]) ):
            x = np.array([bSpline.coefs[0][-1]])
        elif( x < (bSpline.coefs[0][0]) ):
            x = np.array([bSpline.coefs[0][0]])
            #raise ValueError("x is outside bounds of B-spline")
        else:
            x = np.array([x])
        
    (y, dydx) = geometryC.bSplineGeometry(bSpline.knots,bSpline.coefs,x)
    
    return (y, dydx)
    
#==============================================================================
# Calculate volume of axisymmetric nozzle wall using trapezoidal integration
#==============================================================================
def wallVolume(innerWall,thickness):
    
    xVolume = np.linspace(0,innerWall.length,1000)
    volumeIntegrand = np.pi*innerWall.diameter(xVolume)*                     \
      thickness.radius(xVolume) + np.pi*thickness.radius(xVolume)**2
    volume = scipy.integrate.trapz(volumeIntegrand,xVolume)
    
    return volume
