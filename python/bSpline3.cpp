#include <iostream>
#include <math.h>
#include <stdlib.h>

extern "C" {

/* Find first 1-based index where scalar xFind < xVec[i] */
int find(double xFind, double *xVec, int size) {
  int counter = 1;
  for(int ii = 1; ii < size; ii++) {
    if(xFind < xVec[ii]) {break;} // assume xVec is in ascending order
    counter += 1;
  }
  return counter;
}

/* Calculate x, y, dxdu, and dydu given scalar value of parameter u for 3rd 
   degree NURB spline */
void uMap3(double *knots, double *coefs, double u, double *x, double *y,
	   double *dxdu, double *dydu, int k, int c)
{

  // Zero all outputs
  *x = 0.;
  *y = 0.;
  *dxdu = 0.;
  *dydu = 0.;

  // Calculate the highest knot index that gives a value of u below
  // the given u value (1-based index)
  int hh = find(u,knots,k);
  if(hh == k) { // at end of knots vector
    hh = k - 4;
  }
  
  int nn;
  if(hh == 1) {nn = 1;}
  else if(hh == 2) {nn = 2;}
  else if(hh == 3) {nn = 3;}
  else {nn = 4;}
  
  int ii = 0;
  int jj;
  double k1, k2, k3, k4, k5; // knot values
  double Ncurrent, dNducurrent; // basis contributions
  
  // Sum up contributions from each basis
  while(ii > -nn) { // i.e. for each contributing basis
    
    jj = hh + ii - 1;
    
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
        Ncurrent = (u-k1)/(k4-k1)*((u-k1)/(k3-k1)*(k3-u)/(k3-k2) + (k4-u)/
          (k4-k2)*(u-k2)/(k3-k2)) + (k5-u)/(k5-k2)*(u-k2)/(k4-k2)*(u-k2)/
          (k3-k2);
        dNducurrent = (((k1 - u)*(k3 - u))/((k1 - k3)*(k2 - k3)) + ((k2 - u)*
          (k4 - u))/((k2 - k3)*(k2 - k4)))/(k1 - k4) + pow(k2 - u,2)/((k2 - k3)*
          (k2 - k4)*(k2 - k5)) + ((k5 - u)*(2*k2 - 2*u))/((k2 - k3)*(k2 - k4)*
          (k2 - k5)) + (2*(k1 - u)*(k1*k2 - k3*k4 - k1*u - k2*u + k3*u + k4*u))/
          ((k1 - k3)*(k1 - k4)*(k2 - k3)*(k2 - k4));
      }    
    } else if(ii == -2) { // calculate basis N3
      if( fabs(k3-k4) <= 1e-8) {Ncurrent = 0; dNducurrent = 0;}
      else {
        Ncurrent = (u-k1)/(k4-k1)*(k4-u)/(k4-k2)*(k4-u)/(k4-k3) + (k5-u)/
          (k5-k2)*((u-k2)/(k4-k2)*(k4-u)/(k4-k3) + (k5-u)/(k5-k3)*(u-k3)/
          (k4-k3));
        dNducurrent = - (((k2 - u)*(k4 - u))/((k2 - k4)*(k3 - k4)) + ((k3 - u)*
          (k5 - u))/((k3 - k4)*(k3 - k5)))/(k2 - k5) - pow(k4 - u,2)/((k1 - k4)*
          (k2 - k4)*(k3 - k4)) - ((k1 - u)*(2*k4 - 2*u))/((k1 - k4)*(k2 - k4)*
          (k3 - k4)) - (2*(k5 - u)*(k2*k3 - k4*k5 - k2*u - k3*u + k4*u + k5*u))/
          ((k2 - k4)*(k2 - k5)*(k3 - k4)*(k3 - k5));
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

/* u is a scalar value of the B-spline parameteric parameter u, knots is 
   a 1-D array of size k, coefs is a 1-D array of size 2*c, where x
   coordinates are listed first, followed by y coordinates */
void bSplineGeo3(double *knots, double *coefs, double *x, double *y,
		     double *dydx, int nx, int k, int c)
{

  /* Check degree of spline */
  int p = k - c - 1;
  if(p != 3) {
    std::cout << "Only B-splines of degree 3 are implemented"
	      << " (degree " << p << " given)" << std::endl;
    return;
  }

  /* Check knots vector format */
  /* For 3rd degree B-spline that starts at first coefficient and ends at last coefficient,
     first 4 and last 4 knots should be duplicated */
  if(knots[0] != knots[1] || knots[1] != knots[2] || knots[2] != knots[3]) {
    std::cout << "First 4 knots should be the same" << std::endl;
    return;
  }
  if(knots[k-1] != knots[k-2] || knots[k-2] != knots[k-3] || knots[k-3] != knots[k-4]) {
    std::cout << "Last 4 knots should be the same" << std::endl;
    return;
  }
  
  /* Check coefficients format?? */

  double xTemp, yTemp, dxduTemp, dyduTemp;

  // Determine x-value at breaks
  //double xKnot[k];
  double *xKnot;
  xKnot = (double *)malloc(sizeof(double)*k);
  for(int ii = 0; ii < k; ii++) {
    uMap3(knots,coefs,(double)knots[ii],&xTemp,&yTemp,&dxduTemp,&dyduTemp,k,c);
    xKnot[ii] = xTemp;
  }
  // assume same 4 knots at end of knot vector; so finding upper bound on u works:
  xKnot[k-4] += 1e-6;

  double tolerance = 1e-6; // tolerance for Newton solver

  int seg; // tally variable for segment number
  double uLower, uUpper; // lower and upper bounds on u
  double u, uNew; // values of u used in Newton iterations
  double xEst, dxduEst; // estimated values of x and dxdu
  int counter; // used to terminate Newton iterations
  double errorMeasure; // used to record error in estimate of u
  for(int ii = 0; ii < nx; ii++) {

    /* Check that x is within proper range */
    if(x[ii] > coefs[c-1]) {
      std::cout << "x (" << x[ii] << ") not within range specified by coefs vector"
		<< std::endl;
      return;
    }
    if(x[ii] < coefs[0]) {
      std::cout <<  "x (" << x[ii] << ") not within range specified by knots vector"
		<< std::endl;
      return;
    }

    // Determine lower and upper bounds on u
    seg = find(x[ii],xKnot,k);
    uLower = knots[seg - 1];
    uUpper = knots[seg];
      
    // Pick a guess for u (a linear interpolation)
    //u = (x[ii] - xKnot[seg-1])/(xKnot[seg] - xKnot[seg-1])*(uUpper - uLower) + uLower;
    u = (uLower + uUpper)/2;
      
    // Calculate x and dxdu corresponding to u
    uMap3(knots,coefs,u,&xTemp,&yTemp,&dxduTemp,&dyduTemp,k,c);
    xEst = xTemp;
    dxduEst = dxduTemp;
      
    // Perform 1 Newton iteration
    if(dxduEst < tolerance) { uNew = 0.; }
    else { uNew = u - (xEst - x[ii])/dxduEst; }
      
    // Perform remaining Newton iterations
    counter = 0;
    errorMeasure = fabs((uNew-u)/uNew);
    while( errorMeasure > tolerance ) {
      
      u = uNew;
       
      uMap3(knots,coefs,u,&xTemp,&yTemp,&dxduTemp,&dyduTemp,k,c);
      xEst = xTemp;
      dxduEst = dxduTemp;       
        
      if(dxduEst < 1e-12) { uNew = u; }
      else { uNew = u - (xEst - x[ii])/dxduEst; }
        
      counter = counter + 1;
       
      if( counter > 20) { break; }

      errorMeasure = fabs((uNew-u)/uNew);
       
    }
      
    u = uNew;
    uMap3(knots,coefs,u,&xTemp,&yTemp,&dxduTemp,&dyduTemp,k,c);
      
    y[ii] = yTemp;
      
    if( dxduTemp < tolerance ) { dydx[ii] = 0; }
    else { dydx[ii] = dyduTemp/dxduTemp; }      
    
  }

  // Free dynamically allocated memory
  free(xKnot);

  return;

}

} // extern "C"
