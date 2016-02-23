function [  ] = nozzleCFDGmsh( nozzle, xwall, ywall )
	% Victorien Menier, Feb 2016
	
	DatOut=fopen('axinoz.geo','w');

	%---- Define mesh element sizes
	
	sizWal = nozzle.sizWal;
	sizFar = nozzle.sizFar;
	sizSym = nozzle.sizSym;
	
	nbv=length(xwall);
	
	CrdBox=[];
		
	CrdBox(1,1) = 0;       CrdBox(1,2) = 0;
	CrdBox(2,1) = 0.67;    CrdBox(2,2) = 0;
	CrdBox(3,1) = 20.1;    CrdBox(3,2) = 0;
	CrdBox(4,1) = 20.1;    CrdBox(4,2) = 6.03;
	CrdBox(5,1) = -0.67;   CrdBox(5,2) = 3.015;
	CrdBox(6,1) = -0.67;   CrdBox(6,2) = 0.4244;
	CrdBox(7,1) = 0.1548;  CrdBox(7,2) = 0.4244;
	CrdBox(8,1) = 0.67;    CrdBox(8,2) =  0.3151;
	
	
	%--- Add points
	
	%fprintf(DatOut,'Point(1) = {0, 0, 0, %f};\n', sizWal);
	%fprintf(DatOut,'Point(2) = {0.67, 0, 0, %f};\n', sizWal);
	%fprintf(DatOut,'Point(3) = {30, 0, 0, %f};\n', sizSym);
	%fprintf(DatOut,'Point(4) = {30, 9, 0, %f};\n', sizFar);
	%fprintf(DatOut,'Point(5) = {-1, 4.5, 0, %f};\n', sizFar);
	%fprintf(DatOut,'Point(6) = {-1, 0.768, 0, %f};\n', sizWal);
	%fprintf(DatOut,'Point(7) = {-0.52, 0.768, 0, %f};\n', sizWal);
	%fprintf(DatOut,'Point(8) = {0.67, 0.37, 0, %f};\n', sizWal);
	
	%for i=1:8
	%	fprintf(DatOut,'Point(%d) = {%f, %f, 0, %f};\n', i, CrdBox(i,1), CrdBox(i,2), );
	%end
	
	fprintf(DatOut,'Point(1) = {%f, %f, 0, %f};\n',  CrdBox(1,1), CrdBox(1,2), sizWal);
	fprintf(DatOut,'Point(2) = {%f, %f, 0, %f};\n',  CrdBox(2,1), CrdBox(2,2), sizWal);
	fprintf(DatOut,'Point(3) = {%f, %f, 0, %f};\n',  CrdBox(3,1), CrdBox(3,2), sizSym);
	fprintf(DatOut,'Point(4) = {%f, %f, 0, %f};\n',  CrdBox(4,1), CrdBox(4,2), sizFar);
	fprintf(DatOut,'Point(5) = {%f, %f, 0, %f};\n',  CrdBox(5,1), CrdBox(5,2), sizFar);
	fprintf(DatOut,'Point(6) = {%f, %f, 0, %f};\n',  CrdBox(6,1), CrdBox(6,2), sizWal);
	fprintf(DatOut,'Point(7) = {%f, %f, 0, %f};\n',  CrdBox(7,1), CrdBox(7,2), sizWal);
	fprintf(DatOut,'Point(8) = {%f, %f, 0, %f};\n',  CrdBox(8,1), CrdBox(8,2), sizWal);
	
	
	%--- Add B-Spline control points
	
	vid=8;
	for i = nbv:-1:1
		vid = vid+1;
		fprintf(DatOut,'Point(%d) = {%f, %f, 0, %f};\n', vid, xwall(i), ywall(i), sizWal);
	end	

	%--- Add lines
	
	fprintf(DatOut,'Line(1) = {1, 2};\n');
	fprintf(DatOut,'Line(2) = {2, 3};\n');
	fprintf(DatOut,'Line(3) = {3, 4};\n');
	fprintf(DatOut,'Line(4) = {4, 5};\n');
	fprintf(DatOut,'Line(5) = {5, 6};\n');
	fprintf(DatOut,'Line(6) = {6, 7};\n');
	fprintf(DatOut,'Line(7) = {7, 8};\n');
	fprintf(DatOut,'Line(8) = {8, 9};\n');
	
	%--- B-Spline
	
	fprintf(DatOut,'BSpline(9) = { 9');
	eid=9;
	
	for vid = 10:10+nbv-2
		fprintf(DatOut,', %d', vid);
	end
	
	fprintf(DatOut,'};\n');
	
	fprintf(DatOut,'Line(10) = {%d, 1};\n', 8+nbv);
	
	%--- Plane surface
	
	fprintf(DatOut,'Line Loop(14) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};\n');
	fprintf(DatOut,'Plane Surface(14) = {14};\n');
	
	if (strcmp(nozzle.governing,'euler'))
		fprintf(DatOut,'Line(11) = {9, 2};\n');
		fprintf(DatOut,'Line{11} In Surface{14};\n');
	elseif (strcmp(nozzle.governing,'rans'))
		
		% Compute ds
		% Cf http://www.pointwise.com/yplus/
		
		Rex  = nozzle.boundaryCdt.Re;
		rho  = nozzle.boundaryCdt.RhoRef;
		Uinf = nozzle.boundaryCdt.Uref;
		yplus = nozzle.yplus;
		mu   = nozzle.boundaryCdt.MuRef;
		
		Cf   = 0.026/Rex^(1./7.);
		tauw = 0.5*Cf*rho*Uinf*Uinf;
		Ufric = sqrt(tauw/rho);
		ds = yplus*mu/(Ufric*rho);
		
		fprintf(DatOut,'Field[1] = BoundaryLayer;\n');
		fprintf(DatOut,'Field[1].EdgesList = {6,7,8,9};\n');
		fprintf(DatOut,'Field[1].NodesList = {6,%d};\n',vid);
		fprintf(DatOut,'Field[1].hfar = 1;\n');
		fprintf(DatOut,'Field[1].hwall_n = %f;\n', ds);
		fprintf(DatOut,'Field[1].hwall_t = %f;\n', sizWal);
		fprintf(DatOut,'Field[1].ratio = 1.3;\n');
		fprintf(DatOut,'Field[1].thickness = 0.02;\n');
		fprintf(DatOut,'BoundaryLayer Field = 1;\n');
	end
	
	
	fclose(DatOut);
	
end