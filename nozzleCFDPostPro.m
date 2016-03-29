
function [nozzle] = nozzleCFDPostPro(meshSU2, Sol, nozzle, fluid, freestream)
	
	%[NbrVer,dat,namVar] = ReadSU2Sol(SolNam);
	
	NbrVer = Sol.NbrVer;
	dat    = Sol.Dat;
	namVar = Sol.NamVar;
	
	NbrVar = size(namVar,2);
	xextract = nozzle.geometry.length; % x coordinate of the exit
    

	thrust = 0;
    mdot = 0; % mass flow rate
	
	ix     = -1;
	iMach  = -1;
	iTem   = -1;
	iRho   = -1;
	iCons1 = -1;
	iCons2 = -1;
	iCons3 = -1;
	iCons4 = -1;
	iPres  = -1;
	iVid   = -1;
	
	for i=1:NbrVar
		if ( strcmp(namVar(i),'"x"') ) 
			ix = i;
		elseif ( strcmp(namVar(i),'"y"') ) 
			iy = i;
		elseif ( strcmp(namVar(i),'"Mach"') ) 
			iMach = i;
		elseif ( strcmp(namVar(i),'"PointID"') ) 
			iVid = i;
		elseif ( strcmp(namVar(i),'"Temperature"') ) 
			iTem = i;
		elseif ( strcmp(namVar(i),'"Conservative_1"') ) 
			iCons1 = i;
		elseif ( strcmp(namVar(i),'"Conservative_2"') ) 
			iCons2 = i;
		elseif ( strcmp(namVar(i),'"Conservative_3"') ) 
			iCons3 = i;
		elseif ( strcmp(namVar(i),'"Conservative_4"') ) 
			iCons4 = i;
		elseif ( strcmp(namVar(i),'"Pressure"') ) 
			iPres = i;
		end
	end
		
	% --------------------------------------------------------
	% ---    Interpolate solution to structured grid
	% --------------------------------------------------------
	
	Ni=10;
	Nj=10;
	
	x = dat(:,ix);
	y = dat(:,iy);
	v = dat(:,iTem);
	
	xq = repmat(linspace(0,nozzle.geometry.length,Ni),Nj,1);
	yq = zeros(Nj,Ni);
	
	[ A, dAdx, D, nozzle ] = nozzleParameterization( nozzle );
	t = @(x) piecewiseLinearGeometry(x,'t',nozzle.wall); % m, thickness of wall
	
	for i=1:size(xq,2)
		yq(:,i) = linspace(0, 0.5*D(xq(1,i)), Nj);
	end
	
	% --------------------------------------------------------
	% ---    Extract solution along line(s) of interest
	% --------------------------------------------------------
	
	idx = find( dat(:,ix) == xextract & dat(:,iy) < 0.2919 );
	DatLin = dat(idx,:);
	DatLin = sortrows(DatLin,iy);
	
	NbvLin = size(idx,1);
	
	if ( NbvLin < 2 ) 
		error('  ## ERROR nozzleCFDPostPro : No extraction data was found');
	end 
	
	avgDat=zeros(1,size(DatLin,2));
  
  for i=2:NbvLin
  	dy = DatLin(i,iy)-DatLin(i-1,iy);
  	avgDat(1,:) = avgDat(1,:)+dy*DatLin(i,:);
		
		rhoU = DatLin(i,iCons2);
		rho  = DatLin(i,iCons1);
		U    = DatLin(i,iCons2)/DatLin(i,iCons1);
		U0   = freestream.U;
		P    = DatLin(i,iPres);
		P0   = freestream.P;
		
		% --- Compute thrust
		% T = 2PI * Int_{0}^{R} (rho U ( U - U0) + P - Po ) r dr
		thrust = thrust + dy*(rhoU*(U-U0)+P-P0);
        mdot = mdot + dy*rhoU;
  end
		
	hei = DatLin(NbvLin,iy)-DatLin(1,iy);
	avgDat(1,:) = avgDat(1,:)/hei;
	
	% --- Compute hoop stress
	
	P = dat(:,iPres);
	Pq = griddata(x,y,P,xq,yq);
	
	nozzle.hoopStress = prod([Pq(Nj,:) ;D(xq(1,:)')';1./(2*t(xq(1,:)')')]);

	% --- Compute
	%WriteGMFStruc ('presin.mesh', xq, yq, Pq )
	
	% --------------------
	% ---    Outputs      
	% --------------------
	
	nozzle.exit.M = avgDat(1,iMach);
	nozzle.exit.T = avgDat(1,iTem);
	nozzle.exit.U = (avgDat(1,iCons2)+avgDat(1,iCons3))/avgDat(1,iCons1);
	nozzle.exit.P = avgDat(1,iPres);
	nozzle.exit.Pstag = nozzle.exit.P*(1 + (fluid.gam-1)*nozzle.exit.M^2/2)^(fluid.gam/(fluid.gam-1));

  nozzle.exit.Tstag = nozzle.exit.T*(1 + (fluid.gam-1)*nozzle.exit.M^2/2);
  nozzle.netThrust = thrust;
  nozzle.massFlowRate = mdot;
  
	fprintf('\n -- Info: CFD results (averaged values at nozzle exit)\n');
	fprintf('           Mach   = %f\n', nozzle.exit.M );
	fprintf('           Temp   = %f K\n', nozzle.exit.T );
	fprintf('           Vel    = %f m/s\n', nozzle.exit.U );
	fprintf('           Pres   = %f Pa\n', nozzle.exit.P );
	fprintf('           Thrust = %f N\n', nozzle.netThrust );
	
end