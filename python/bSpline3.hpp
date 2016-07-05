#ifndef BSPLINE3_HPP
#define BSPLINE3_HPP

extern "C" {

/* Forward search through array for first element greater than a given xFind */
int find(double xFind, double *xVec, int size);

/* 3rd degree B-spline mapping from parametric u to coordinates x, y, and 
   derivatives */
void uMap3(double *knots, double *coefs, double u, double *x, double *y,
	   double *dxdu, double *dydu, int k, int c);

/* Given x and B-spline parameters, calculate y and dydx for 2-D 3rd degree 
   B-spline */
void bSplineGeo3(double *knots, double *coefs, double *x, double *y,
		 double *dydx, int nx, int k, int c);

}

#endif /* BSPLINE3_HPP */
