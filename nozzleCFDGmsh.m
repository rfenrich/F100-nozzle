function [  ] = nozzleCFDGmsh( nozzle, xwall, ywall )
	% Victorien Menier, Feb 2016
	
	DatOut=fopen('axinoz.geo','w');

	%---- Define mesh element sizes
	
	sizWal = nozzle.sizWal;
	sizFar = nozzle.sizFar;
	sizSym = nozzle.sizSym;
	
	sizWal = 0.3;
	sizFar = 0.3;
	sizSym = 0.3;
	
	nbv=length(xwall);
	
	CrdBox=[];
		
	%CrdBox(1,1) = 0;       CrdBox(1,2) = 0;
	%CrdBox(2,1) = 0.67;    CrdBox(2,2) = 0;
	%CrdBox(3,1) = 20.1;    CrdBox(3,2) = 0;
	%CrdBox(4,1) = 20.1;    CrdBox(4,2) = 6.03;
	%CrdBox(5,1) = -0.67;   CrdBox(5,2) = 3.015;
	%CrdBox(6,1) = -0.67;   CrdBox(6,2) = 0.4244;
	%CrdBox(7,1) = 0.1548;  CrdBox(7,2) = 0.4244;
	%CrdBox(8,1) = 0.67;    CrdBox(8,2) =  0.3151;
	
	% --- Updated domain : shorter and higher (to prevent unwanted physics) computational domain
	CrdBox(1,1) = 0;                         CrdBox(1,2) = 0;
	CrdBox(2,1) = nozzle.geometry.length;    CrdBox(2,2) = 0;
	CrdBox(3,1) = 1.5;                         CrdBox(3,2) = 0;
	CrdBox(4,1) = 1.5;                         CrdBox(4,2) = 2.5;
	CrdBox(5,1) = -0.67;                     CrdBox(5,2) = 2.5;
	CrdBox(6,1) = -0.67;                     CrdBox(6,2) = 0.4244;
	CrdBox(7,1) = 0.1548;                    CrdBox(7,2) = 0.4244;
	CrdBox(8,1) = nozzle.geometry.length;    %CrdBox(8,2) =  0.3151;
	CrdBox(8,2) = ywall(end)+0.012;
	
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
	
	fprintf(DatOut,'Point(9)  = {%f, %f, 0, %f};\n',  CrdBox(3,1), CrdBox(8,2), sizWal);
	fprintf(DatOut,'Point(10) = {%f, %f, 0, %f};\n',  CrdBox(3,1), CrdBox(7,2)+0.25*CrdBox(8,2), sizWal);
	fprintf(DatOut,'Point(11) = {%f, %f, 0, %f};\n',  CrdBox(6,1), CrdBox(7,2)+0.25*CrdBox(8,2), sizWal);
	
	%--- Add B-Spline control points
	
	vid=11;
	for i = nbv:-1:1
		vid = vid+1;
		fprintf(DatOut,'Point(%d) = {%f, %f, 0, %f};\n', vid, xwall(i), ywall(i), sizWal);
	end	

	%--- Add lines
	
	fprintf(DatOut,'Line(1)  = {1, 2};\n');
	fprintf(DatOut,'Line(2)  = {2, 3};\n');
	fprintf(DatOut,'Line(3)  = {3, 9};\n');
	fprintf(DatOut,'Line(4)  = {9, 10};\n');
	fprintf(DatOut,'Line(5)  = {10, 4};\n');
	fprintf(DatOut,'Line(6)  = {4, 5};\n');
	fprintf(DatOut,'Line(7)  = {5, 11};\n');
	fprintf(DatOut,'Line(8)  = {11, 6};\n');
	fprintf(DatOut,'Line(9)  = {6, 7};\n');
	fprintf(DatOut,'Line(10) = {7, 8};\n');
	fprintf(DatOut,'Line(11) = {8, 12};\n');
				
	
	%--- B-Spline
	
	fprintf(DatOut,'BSpline(12) = { 12');
	eid=12;
	
	for vid = eid:eid+nbv-1
		fprintf(DatOut,', %d', vid);
	end
	
	fprintf(DatOut,'};\n');
	
	fprintf(DatOut,'Line(13) = {%d, 1};\n', vid);
	
	%--- Plane surface
	
	fprintf(DatOut,'Line Loop(14) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};\n');
	fprintf(DatOut,'Plane Surface(14) = {14};\n');
	
	if (strcmp(nozzle.governing,'euler'))
		%fprintf(DatOut,'Line(11) = {9, 2};\n');
		%fprintf(DatOut,'Line{11} In Surface{14};\n');
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
	
	
	%Mesh.Algorithm = 5;
	
	NbrFld = 6;
	fields = zeros(NbrFld,5);
	%volume boundingbox  -100 100  -10 0.8  -100 100 h 0.05
	%volume boundingbox  -100 100  -10 0.5  -100 100 h 0.01
	%volume boundingbox  -100 100  -10 0.3048  -100 100 h 0.005
	%volume boundingbox  -100 0.3  -10 0.37  -100 100 h 0.005
	%volume boundingbox  -100 100  -10 0.8  -100 100 h 0.07
	
	hl1 = 0.1;
	hl2 = 0.07;
	hl3 = 0.05;
	hl4 = 0.005;
	hl5 = 0.009;
	
	hl1 = hl1;
	hl2 = hl2;
	hl3 = 1.2*hl3;
	hl4 = 1.2*hl4;
	hl5 = 1.2*hl5;
	
	%volume boundingbox  -100 100  -10 0.8  -100 100 h 0.05
	%volume boundingbox  -100 100  -10 0.5  -100 100 h 0.01
	%volume boundingbox  -100 0.7  -10 0.5  -100 100 h 0.005
	%volume boundingbox  -100 100  -10 0.3048  -100 100 h 0.009
	%volume boundingbox  -100 0.3  -10 0.37  -100 100 h 0.005
	%volume boundingbox  -100 100  -10 0.8  -100 100 h 0.07
	
	% xmin, xmax, ymin, ymax, size
	fields(1,:) = [-100, 100, -10, 0.8,    hl3];
	fields(2,:) = [-100, 100, -10, 0.5,    hl1];
	fields(3,:) = [-100, 100, -10, 0.5,    hl5];
	fields(4,:) = [-100, 0.3, -10, 0.37,   hl4];
	fields(5,:) = [-100, 100, -10, 0.8,    hl2];
	fields(6,:) = [-100, 0.7, -10, 0.5,    hl4];

	
	for i=1:NbrFld
		fprintf(DatOut,'Field[%d] = Box;\n'    ,i);
		fprintf(DatOut,'Field[%d].VIn = %f;\n' ,i,fields(i,5));
		fprintf(DatOut,'Field[%d].VOut = %f;\n',i,sizFar);
		fprintf(DatOut,'Field[%d].XMin = %f;\n',i,fields(i,1));
		fprintf(DatOut,'Field[%d].XMax = %f;\n',i,fields(i,2));
		fprintf(DatOut,'Field[%d].YMin = %f;\n',i,fields(i,3));
		fprintf(DatOut,'Field[%d].YMax = %f;\n',i,fields(i,4));
	end
	

	fprintf(DatOut,'Field[%d] = Min;\n', NbrFld+1);
	fprintf(DatOut,'Field[%d].FieldsList = { 1', NbrFld+1);	
	for i=2:NbrFld
		fprintf(DatOut,', %d ', i);
	end
	fprintf(DatOut,'};\n');

	fprintf(DatOut,'Background Field = %d;\n', NbrFld+1);
	
	
	
	
	
	fclose(DatOut);
	
end