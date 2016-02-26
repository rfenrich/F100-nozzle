#include "mex.h"
#include <math.h>
#include <vector>

/* Given u calculate y for 2nd degree NURB spline */
double u2y(double *coefs, double *knots, mwSize i, mwSize n, double u) {
  return -(((u - knots[i])*(u - knots[2 + i]))/((knots[i] - knots[2 + i])*(knots[1 + i] - knots[2 + i])) + ((u - knots[1 + i])*(u - knots[3 + i]))/((knots[1 + i] - knots[2 + i])*(knots[1 + i] - knots[3 + i])))*coefs[n + i] + (pow(u - knots[2 + i],2)*coefs[n + i - 1])/((knots[i] - knots[2 + i])*(knots[1 + i] - knots[2 + i])) + (pow(u - knots[1 + i],2)*coefs[n + i + 1])/((knots[1 + i] - knots[2 + i])*(knots[1 + i] - knots[3 + i]));
}

/* Given u calculate x for 2nd degree NURB spline */
double u2x(double *coefs, double *knots, mwSize i, mwSize n, double u) {
  return -(((u - knots[i])*(u - knots[2 + i]))/((knots[i] - knots[2 + i])*(knots[1 + i] - knots[2 + i])) + ((u - knots[1 + i])*(u - knots[3 + i]))/((knots[1 + i] - knots[2 + i])*(knots[1 + i] - knots[3 + i])))*coefs[i] + (pow(u - knots[2 + i],2)*coefs[i - 1])/((knots[i] - knots[2 + i])*(knots[1 + i] - knots[2 + i])) + (pow(u - knots[1 + i],2)*coefs[i + 1])/((knots[1 + i] - knots[2 + i])*(knots[1 + i] - knots[3 + i]));
}

/* Given u calculate derivative of g = f(u) w.r.t. u for 2nd degree NURB spline */
double dgdu(double *coefs, double *knots, mwSize i, mwSize n, double u) {
  return -(2*(u*knots[i]*coefs[i] - u*knots[1 + i]*coefs[i - 1] - u*knots[i]*coefs[i + 1] + u*knots[1 + i]*coefs[i] - u*knots[2 + i]*coefs[i] + u*knots[3 + i]*coefs[i - 1] + u*knots[2 + i]*coefs[i + 1] - u*knots[3 + i]*coefs[i] - knots[i]*knots[1 + i]*coefs[i] + knots[i]*knots[1 + i]*coefs[i + 1] + knots[1 + i]*knots[2 + i]*coefs[i - 1] - knots[1 + i]*knots[2 + i]*coefs[i + 1] - knots[2 + i]*knots[3 + i]*coefs[i - 1] + knots[2 + i]*knots[3 + i]*coefs[i]))/((knots[i] - knots[2 + i])*(knots[1 + i] - knots[2 + i])*(knots[1 + i] - knots[3 + i]));
}

/* Given u calculate g = f(u) - x w.r.t. u for 2nd degree NURB spline */
double g(double *coefs, double *knots, mwSize i, mwSize n, double u, double r) {
  return u2x(coefs,knots,i,n,u) - r;
 }
 
/* Given u calculate derivative of y w.r.t. u for 2nd degree NURB spline */
double dydu(double *coefs, double *knots, mwSize i, mwSize n, double u) {
  return -(2*(u*knots[i]*coefs[n + i] - u*knots[1 + i]*coefs[n + i - 1] - u*knots[i]*coefs[n + i + 1] + u*knots[1 + i]*coefs[n + i] - u*knots[2 + i]*coefs[n + i] + u*knots[3 + i]*coefs[n + i - 1] + u*knots[2 + i]*coefs[n + i + 1] - u*knots[3 + i]*coefs[n + i] - knots[i]*knots[1 + i]*coefs[n + i] + knots[i]*knots[1 + i]*coefs[n + i + 1] + knots[1 + i]*knots[2 + i]*coefs[n + i - 1] - knots[1 + i]*knots[2 + i]*coefs[n + i + 1] - knots[2 + i]*knots[3 + i]*coefs[n + i - 1] + knots[2 + i]*knots[3 + i]*coefs[n + i]))/((knots[i] - knots[2 + i])*(knots[1 + i] - knots[2 + i])*(knots[1 + i] - knots[3 + i]));
}

/* Given u calculate derivative of u w.r.t. x for 2nd degree NURB spline */
double dudx(double *coefs, double *knots, mwSize i, mwSize n, double u) {
  return 1/(-(2*(u*knots[i]*coefs[i] - u*knots[1 + i]*coefs[i - 1] - u*knots[i]*coefs[i + 1] + u*knots[1 + i]*coefs[i] - u*knots[2 + i]*coefs[i] + u*knots[3 + i]*coefs[i - 1] + u*knots[2 + i]*coefs[i + 1] - u*knots[3 + i]*coefs[i] - knots[i]*knots[1 + i]*coefs[i] + knots[i]*knots[1 + i]*coefs[i + 1] + knots[1 + i]*knots[2 + i]*coefs[i - 1] - knots[1 + i]*knots[2 + i]*coefs[i + 1] - knots[2 + i]*knots[3 + i]*coefs[i - 1] + knots[2 + i]*knots[3 + i]*coefs[i]))/((knots[i] - knots[2 + i])*(knots[1 + i] - knots[2 + i])*(knots[1 + i] - knots[3 + i])));
}

/* Find first 1-based index where scalar xFind < xVec[i] */
mwSize find(double xFind, double *xVec, mwSize sz) {
  mwSize counter = 1;
  for(mwSize ii = 1; ii < sz; ii++) {
    if(xFind < xVec[ii]) {break;} // assume xVec is in ascending order
    counter += 1;
  }
  return counter;
}

/* Calculate x, y, dxdu, and dydu given u for 3rd degree NURB spline */
/* Currently implemented only for scalar inputs of u */
void uMap3(double *knots, double *coefs, double u, double *x, double *y, double* dxdu, double *dydu, mwSize k, mwSize c) {

  // Zero all outputs
  *x = 0.;
  *y = 0.;
  *dxdu = 0.;
  *dydu = 0.;
    
  // Calculate the highest knot index that gives a value of u below
  // the given value (1-based index)
  mwSize hh = find(u,knots,k);
  if(hh == k) { // at end of knots vector
    hh = k - 4;
  }
  //mexPrintf("hh: %i\n",hh);
  //mexPrintf("u: %f\n",u);
  
  mwSize nn;
  if(hh == 1) {nn = 1;}
  else if(hh == 2) {nn = 2;}
  else if(hh == 3) {nn = 3;}
  else {nn = 4;}
  
  mwSize ii = 0;
  mwSize jj;
  double k1, k2, k3, k4, k5; // knot values
  double Ncurrent, dNducurrent; // basis contributions
  
  // Sum up contributions from each basis
  while(ii > -nn) { // i.e. for each contributing basis
    
    jj = hh + ii - 1; // subtracted 2 for C++ compatibility
    //mexPrintf("jj = %i\n",jj);
    
    // Redefine k1 through k5 here
    k1 = knots[jj];
    k2 = knots[jj+1];
    k3 = knots[jj+2];
    k4 = knots[jj+3];
    k5 = knots[jj+4];
    
    if(ii == 0) { // calculate basis N1
      if( fabs(k1-k2) <= 1e-8) {Ncurrent = 0; dNducurrent = 0;}
      else {
        Ncurrent = (u-k1)/(k4-k1)*(u-k1)/(k3-k1)*(u-k1)/(k2-k1);
        dNducurrent = -(3*pow(k1 - u,2))/((k1 - k2)*(k1 - k3)*(k1 - k4));
      }
    } else if(ii == -1) { // calculate basis N2
      if( fabs(k2-k3) <= 1e-8) {Ncurrent = 0; dNducurrent = 0;}
      else {
        Ncurrent = (u-k1)/(k4-k1)*((u-k1)/(k3-k1)*(k3-u)/(k3-k2) + (k4-u)/(k4-k2)*(u-k2)/(k3-k2)) + (k5-u)/(k5-k2)*(u-k2)/(k4-k2)*(u-k2)/(k3-k2);
        dNducurrent = (((k1 - u)*(k3 - u))/((k1 - k3)*(k2 - k3)) + ((k2 - u)*(k4 - u))/((k2 - k3)*(k2 - k4)))/(k1 - k4) + pow(k2 - u,2)/((k2 - k3)*(k2 - k4)*(k2 - k5)) + ((k5 - u)*(2*k2 - 2*u))/((k2 - k3)*(k2 - k4)*(k2 - k5)) + (2*(k1 - u)*(k1*k2 - k3*k4 - k1*u - k2*u + k3*u + k4*u))/((k1 - k3)*(k1 - k4)*(k2 - k3)*(k2 - k4));
      }    
    } else if(ii == -2) { // calculate basis N3
      if( fabs(k3-k4) <= 1e-8) {Ncurrent = 0; dNducurrent = 0;}
      else {
        Ncurrent = (u-k1)/(k4-k1)*(k4-u)/(k4-k2)*(k4-u)/(k4-k3) + (k5-u)/(k5-k2)*((u-k2)/(k4-k2)*(k4-u)/(k4-k3) + (k5-u)/(k5-k3)*(u-k3)/(k4-k3));
        dNducurrent = - (((k2 - u)*(k4 - u))/((k2 - k4)*(k3 - k4)) + ((k3 - u)*(k5 - u))/((k3 - k4)*(k3 - k5)))/(k2 - k5) - pow(k4 - u,2)/((k1 - k4)*(k2 - k4)*(k3 - k4)) - ((k1 - u)*(2*k4 - 2*u))/((k1 - k4)*(k2 - k4)*(k3 - k4)) - (2*(k5 - u)*(k2*k3 - k4*k5 - k2*u - k3*u + k4*u + k5*u))/((k2 - k4)*(k2 - k5)*(k3 - k4)*(k3 - k5));
      }       
    } else { // calculate basis N4
      if( fabs(k4-k5) <= 1e-8) {Ncurrent = 0; dNducurrent = 0;}
      else {
        Ncurrent = (k5-u)/(k5-k2)*(k5-u)/(k5-k3)*(k5-u)/(k5-k4);  
        dNducurrent = (3*pow(k5 - u,2))/((k2 - k5)*(k3 - k5)*(k4 - k5));
      }     
    }
    
    *x += coefs[jj]*Ncurrent;
    *y += coefs[c + jj]*Ncurrent;
    *dxdu += coefs[jj]*dNducurrent;
    *dydu += coefs[c + jj]*dNducurrent;
    
    ii -= 1;
    
    if(ii < -4) {break;}
    
  }  
  
  return;
}

/* C++ equivalent of Matlab BsplineGeometry.m function, although it will not calculate 
   a new throat location, or output thickness 't' or y-coordinate 'y'. Given a scalar or
   column vector of x-positions in the vector x, this function calculates the area, 
   change in area dAdx, and diameter at each x-location, using the given knots vector and 
   coefficients matrix for a 2nd degree NURB spline. The output is managed by an output 
   flag, where (1 corresponds to area A, 2 corresponds to dAdx, 3 corresponds to diameter 
   D, otherwise all three are output). Requirements: The knot vector is assumed to be an 
   ascending vector with elements spaced by either 0 or 1, and with 3 repeated elements at 
   the beginning and end. For k total knots, there must be k-3 coefficients to achieve a 
   2nd degree NURB spline. The repeated elements in the knot vector ensure the NURB spline
   begins at the first coefficient and ends at the last coefficient. */
void mexFunction(mwSize nlhs, mxArray *plhs[],  /* Output variables */
				 mwSize nrhs, const mxArray *prhs[])    /* Input variables */
{	

  /* First check inputs */

  if(nrhs != 4) {
    mexErrMsgTxt("Four inputs required: x (n x 1 vector), "
                 "outputResult (int, where 1 outputs A, 2 outputs dAdx, 3 outputs D), "
                 "knots (k x 1 vector), coefs (c x 2 matrix).");
  }
  
  /* Check that x is a vector of doubles */
  if( !mxIsDouble(prhs[0]) || 
    mxIsComplex(prhs[0])) {
    mexErrMsgTxt("x must be type double.");
  }
  /* check that x is a column vector */
  if(mxGetN(prhs[0]) != 1) {
    mexErrMsgTxt("x must be a column vector.");
  }
  
  /* Check that outputResult is a scalar */
  if( !mxIsDouble(prhs[1]) || 
     mxIsComplex(prhs[1]) ||
     mxGetNumberOfElements(prhs[1]) != 1 ) {
    mexErrMsgTxt("outputResult must be a scalar.");
  }
  
  /* check that knots is a row vector */
  if(mxGetN(prhs[2]) != 1) {
    mexErrMsgTxt("knots must be a column vector.");
  }
  
  /* Check that coefs is a vector of doubles */
  if( !mxIsDouble(prhs[3]) || 
    mxIsComplex(prhs[3])) {
    mexErrMsgTxt("coefs must be type double.");
  }
  /* check that coefs has only two columns */
  if(mxGetN(prhs[3]) != 2) {
    mexErrMsgTxt("coefs must be a matrix with x-coordinates "
                 "in the first column and y-coordinates in "
                 "the second.");
  }
  
  /* Now assign inputs and other necessary variables */
  
  const double PI = 3.141592653589793;
	
	// Inputs include: x, outputResult, knots, coefs
	double *x = mxGetPr(prhs[0]);
	mwSize outputResult = (mwSize)mxGetScalar(prhs[1]);
	double *knots = mxGetPr(prhs[2]);
	double *coefs = mxGetPr(prhs[3]);

  mwSize nx, c, k, p;
	nx = (mwSize)mxGetNumberOfElements(prhs[0]); // number of x
	c = (mwSize)mxGetNumberOfElements(prhs[3])/2; // number of control points
	k = (mwSize)mxGetNumberOfElements(prhs[2]); // number of knots
	p = k - c - 1; // spline degree	
	
	// Outputs include: A, dAdx, D
  double *A, *dAdx, *D, *xThroat, *yThroat;
  double *xdebug, *Ddebug;
  
	/* Prepare output data */
  switch(outputResult) {
  case 1: // return A
    plhs[0] = mxCreateDoubleMatrix(nx,1,mxREAL);
    A = mxGetPr(plhs[0]);
    break;
  case 2: // return dAdx
    plhs[0] = mxCreateDoubleMatrix(nx,1,mxREAL);
    dAdx = mxGetPr(plhs[0]);
    break;
  case 3: // return D
    plhs[0] = mxCreateDoubleMatrix(nx,1,mxREAL);
    D = mxGetPr(plhs[0]);
    break;
  case 4: // return xThroat and yThroat
    mexErrMsgTxt("Calculations for throat currently not enabled in MEX file.");
    break;
  case 8: // for debugging
    nx = (mwSize)1000;
    plhs[0] = mxCreateDoubleMatrix(nx,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nx,1,mxREAL);
    xdebug = mxGetPr(plhs[0]);
    Ddebug = mxGetPr(plhs[1]);
    break;
  default:
    plhs[0] = mxCreateDoubleMatrix(nx,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nx,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nx,1,mxREAL);
    A = mxGetPr(plhs[0]);
    dAdx = mxGetPr(plhs[1]);
    D = mxGetPr(plhs[2]);
  }
	
  /* Check that spline degree = 2 or 3, since calculations are only valid
     for this degree, and for a knot vector spaced only with 0s or 1s */
  if(p == 2) { 
    if(knots[0] != knots[1] && knots[1] != knots[2] && knots[2] != knots[3]) {
      mexErrMsgTxt("2nd degree B-spline implementation requires "
                   "first 3 knots to be identical.");
    }
    if(knots[k] != knots[k-1] && knots[k-1] != knots[k-2] && knots[k-2] != knots[k-3]) {
      mexErrMsgTxt("2nd degree B-spline implementation requires "
                   "last 3 knots to be identical.");
    }    
  } else if (p == 3) {
    if(knots[0] != knots[1] && knots[1] != knots[2] && knots[2] != knots[3] && knots[3] != knots[4]) {
      mexErrMsgTxt("3rd degree B-spline implementation requires "
                   "first 4 knots to be identical.");
    }
    if(knots[k] != knots[k-1] && knots[k-1] != knots[k-2] && knots[k-2] != knots[k-3] && knots[k-3] != knots[k-4]) {
      mexErrMsgTxt("3rd degree B-spline implementation requires "
                   "last 4 knots to be identical.");
    }
    } else {
		mexErrMsgTxt("Calculations are only for a spline of degree 2 or 3.");
    }
	
	double *xKnot;
	double *y, *dydx;
  y = (double *)mxMalloc(8*nx);
  dydx = (double *)mxMalloc(8*nx);
	
	if(p == 2) { // 2nd degree spline
	
    // Determine x value at breaks
    xKnot = (double *)mxMalloc(8*(k-4));
    mwSize seg;
    for(mwSize ii = 0; ii < k - 4; ii++) { // assumes 3 repeated knots at start and end
      // determine segment number
      if(ii < k - 5) {seg = ii + 1;} 
      else {seg = k - 5;}
      xKnot[ii] = u2x(coefs,knots,seg,c,knots[ii+2]);
    }
    xKnot[k-5] += 1e-6; // so following algorithm works
    
    /* Calculate u given x as well as the corresponding y and dydx */
    double uGuess, u;
    for(mwSize ii = 0; ii < nx; ii++) {
      
      // Determine segment number for a given x
      seg = find(x[ii],xKnot,k-5);
          
      // Solve for u using Newton method
      uGuess = (double)seg - 0.5; // assuming knots spaced by 1
      double ftemp1, ftemp2;
      for(mwSize jj = 0; jj < 15; jj++) {
        ftemp1 = g(coefs,knots,seg,c,uGuess,x[ii]);
        ftemp2 = dgdu(coefs,knots,seg,c,uGuess);
        u = uGuess - ftemp1/ftemp2;
        if( fabs((u-uGuess)/uGuess) <= 1e-6) {break;}
        uGuess = u;
      }
      
      // Calculate y corresponding to given x
      y[ii] = u2y(coefs,knots,seg,c,u);
    
      // Calculate dydx
      dydx[ii] = dydu(coefs,knots,seg,c,u)*dudx(coefs,knots,seg,c,u);  
    
    }
  
  } else { // 3rd degree spline
  
    double xTemp, yTemp, dxduTemp, dyduTemp;
  
    // Determine x value at breaks
    xKnot = (double *)mxMalloc(8*k);
    
    for(mwSize ii = 0; ii < k; ii++) {
      uMap3(knots,coefs,(double)knots[ii],&xTemp,&yTemp,&dxduTemp,&dyduTemp,k,c);
      xKnot[ii] = xTemp;
    }
    xKnot[k-4] += 1e-6; // so finding upper bound on u works (assumes same 4 knots at end of knot vector) 
    
    /*for(mwSize ii = 0; ii < k; ii++) {
      mexPrintf("xKnot[%i] = %f\n",ii,xKnot[ii]);
    }*/
    
    double tolerance = 1e-6; // tolerance for Newton solver
    
    mwSize seg; // tally variable for segment number
    double uLower, uUpper; // lower and upper bounds on u
    double u, uNew; // values of u used in Newton iterations
    double xEst, dxduEst; // estimated values of x and dxdu
    mwSize counter; // used to terminate Newton iterations
    double errorMeasure; // used to record error in estimate of u
    for(mwSize ii = 0; ii < nx; ii++) {
    
      // Determine lower and upper bounds on u
      seg = find(x[ii],xKnot,k);
      uLower = knots[seg - 1];
      uUpper = knots[seg];
      //mexPrintf("x = %f\n",x[ii]);
      //mexPrintf("xLower = %f\n",xKnot[seg-1]);
      //mexPrintf("xUpper = %f\n",xKnot[seg]);
      
      // Pick a guess for u (a linear interpolation)
      //u = (x[ii] - xKnot[seg-1])/(xKnot[seg] - xKnot[seg-1])*(uUpper - uLower) + uLower;
      u = (uLower + uUpper)/2;
      //mexPrintf("u: %f\n",u);
      
      // Calculate x and dxdu corresponding to u
      uMap3(knots,coefs,u,&xTemp,&yTemp,&dxduTemp,&dyduTemp,k,c);
      xEst = xTemp;
      dxduEst = dxduTemp;
      //mexPrintf("xEst: %f\n",xEst);
      //mexPrintf("dxduEst: %f\n",dxduEst);
      
      // Perform 1 Newton iteration
      if(dxduEst < tolerance) { uNew = 0.; }
      else { uNew = u - (xEst - x[ii])/dxduEst; }
      //mexPrintf("uNew: %f\n",uNew);
      //mexPrintf("u: %f\n",u);
      
      // Perform remaining Newton iterations
      counter = 0;
      errorMeasure = fabs((uNew-u)/uNew);
      //mexPrintf("errorMeasure: %f\n",errorMeasure);
      while( errorMeasure > tolerance ) {
      
        u = uNew;
        //mexPrintf("u: %f\n",u);
        
        uMap3(knots,coefs,u,&xTemp,&yTemp,&dxduTemp,&dyduTemp,k,c);
        xEst = xTemp;
        dxduEst = dxduTemp;
        //mexPrintf("xEst: %f\n",xEst);
        //mexPrintf("dxduEst: %f\n",dxduEst);        
        
        if(dxduEst < 1e-12) { uNew = u; }
        else { uNew = u - (xEst - x[ii])/dxduEst; }
        
        counter = counter + 1;
        
        if( counter > 20) { break; }

	errorMeasure = fabs((uNew-u)/uNew);
       
      }
      
      u = uNew;
      //mexPrintf("u: %f\n",u);
      uMap3(knots,coefs,u,&xTemp,&yTemp,&dxduTemp,&dyduTemp,k,c);
      //mexPrintf("u: %f\n",u);
      //mexPrintf("x: %f\n",xTemp);
      //mexPrintf("y: %f\n",yTemp);
      //mexPrintf("dxdu: %f\n",dxduTemp);
      //mexPrintf("dydu: %f\n",dyduTemp);
      
      y[ii] = yTemp;
      
      if( dxduTemp < tolerance ) { dydx[ii] = 0; }
      else { dydx[ii] = dyduTemp/dxduTemp; }      
    
    }
  
  }
  
  // Assign data to outputs
  switch(outputResult) {
  case 1: // return A
    for(mwSize ii = 0; ii < nx; ii++) {
      A[ii] = PI*pow(y[ii],2);
    }
    break;
  case 2: // return dAdx
    for(mwSize ii = 0; ii < nx; ii++) {
      dAdx[ii] = 2*PI*y[ii]*dydx[ii];
    }
    break;
  case 3: // return D
    for(mwSize ii = 0; ii < nx; ii++) {
      D[ii] = 2*y[ii];
    }
    break;
  case 4: // return xThroat and yThroat
    mexErrMsgTxt("Calculations for throat currently not enabled in MEX file.");
    break;
  case 8: // return x and D
    double uTemp, xTemp, yTemp, dxduTemp, dyduTemp;
    uTemp = 0.;
    double du;
    du = (double)knots[k-1]/(nx-1);
    for(mwSize ii = 0; ii < nx; ii++) {
      uMap3(knots,coefs,uTemp,&xTemp,&yTemp,&dxduTemp,&dyduTemp,k,c);
      xdebug[ii] = xTemp;
      Ddebug[ii] = yTemp*2;
      uTemp += du;
    }
    break;
  default:
    for(mwSize ii = 0; ii < nx; ii++) {
      A[ii] = PI*pow(y[ii],2);
      dAdx[ii] = 2*PI*y[ii]*dydx[ii];
      D[ii] = 2*y[ii];
    }
  }
	
	mxFree(xKnot);
	mxFree(y);
	mxFree(dydx);
	
	return;
}
