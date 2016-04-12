
function [nozzle] = nozzleCFDPostPro(meshSU2, Sol, nozzle, fluid, freestream)
	
	%[NbrVer,dat,namVar] = ReadSU2Sol(SolNam);
	
	NbrVer = Sol.NbrVer;
	dat    = Sol.Dat;
	namVar = Sol.NamVar;
	
	NbrVar = size(namVar,2);
	xextract = nozzle.geometry.length; % x coordinate of the exit
    

	thrust = 0;
    mdot = 0; % mass flow rate
    PstagExitAverage = 0;
    TstagExitAverage = 0;
	
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
	
	Ni=200;
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
	
	%idx = find( dat(:,ix) == xextract & dat(:,iy) < 0.2919 );
    idx = find( dat(:,ix) == xextract & dat(:,iy) < nozzle.geometry.D(end)/2 + 1e-6 );
	DatLin = dat(idx,:);
	DatLin = sortrows(DatLin,iy);
	
	NbvLin = size(idx,1);
	
	if ( NbvLin < 2 ) 
		error('  ## ERROR nozzleCFDPostPro : No extraction data was found');
	end 
	
	avgDat=zeros(1,size(DatLin,2));
  
    yTotal = max(DatLin(:,iy));
  for i=2:NbvLin
  	dy = DatLin(i,iy)-DatLin(i-1,iy);
  	avgDat(1,:) = avgDat(1,:)+dy*DatLin(i,:);
    
        rhoU = DatLin(i,iCons2);
		rho  = DatLin(i,iCons1);
		U    = DatLin(i,iCons2)/DatLin(i,iCons1);
		U0   = freestream.U;
		P    = DatLin(i,iPres);
		P0   = freestream.P;
        M    = DatLin(i,iMach);
        T    = DatLin(i,iTem);
           
        % --- Compute stagnation properties at exit
        PstagExitAverage = PstagExitAverage + (dy/yTotal)*P*(1 + (fluid.gam-1)*M^2/2)^(fluid.gam/(fluid.gam-1));
        TstagExitAverage = TstagExitAverage + (dy/yTotal)*T*(1 + (fluid.gam-1)*M^2/2); 
		
		% --- Compute thrust
		% T = 2PI * Int_{0}^{R} (rho U ( U - U0) + P - Po ) r dr
		thrust = thrust + dy*(rhoU*(U-U0)+P-P0);
        mdot = mdot + dy*rhoU;
  end
		
	hei = DatLin(NbvLin,iy)-DatLin(1,iy);
	avgDat(1,:) = avgDat(1,:)/hei;
	
	% --- Compute hoop stress
	
    % --- Extract pressure, temp, density, and velocity onto a grid
	P = dat(:,iPres);
	Pq = griddata(x,y,P,xq,yq);
    T = dat(:,iTem);
    Tq = griddata(x,y,T,xq,yq);
    rho = dat(:,iCons1);rho = dat(:,iCons1);
    rhoq = griddata(x,y,rho,xq,yq);
    U = dat(:,iCons2)./dat(:,iCons1);
    Uq = griddata(x,y,U,xq,yq);
    
	% --- Compute
	%WriteGMFStruc ('presin.mesh', xq, yq, Pq )
	
	% --------------------
	% ---    Outputs      
	% --------------------
    
    % --- Calculate geometry
    nozzle.xPosition = xq(1,:);
    nozzle.geometry.A = A(nozzle.xPosition');
    nozzle.geometry.dAdx = dAdx(nozzle.xPosition');
    nozzle.geometry.D = D(nozzle.xPosition');
    nozzle.wall.t = t(nozzle.xPosition);
    nozzle.geometry.maxSlope = max(nozzle.geometry.dAdx./pi./nozzle.geometry.D);
    nozzle.geometry.minSlope = min(nozzle.geometry.dAdx./pi./nozzle.geometry.D);
    
    % --- Assign flow data to nozzle
    nozzle.flow.T = Tq;
    nozzle.flow.P = Pq;
    nozzle.flow.density = rhoq;
    nozzle.flow.U = Uq;
    nozzle.flow.M = Uq./sqrt(fluid.gam*fluid.R*Tq);
    nozzle.flow.Tstag = nozzle.flow.T.*(1 + (fluid.gam-1).*nozzle.flow.M.^2/2); 
    nozzle.flow.Pstag = nozzle.flow.P.*(1 + (fluid.gam-1).*nozzle.flow.M.^2/2).^(fluid.gam/(fluid.gam-1));
    
    % Sutherland's law for dynamic viscosity:
    dynamicViscosity = @(T) 1.716e-5*(T/273.15).^1.5*(273.15 + 110.4)./(T + 110.4); % kg/m*s
    nozzle.flow.Re = rhoq.*Uq.*repmat(nozzle.geometry.D,1,Nj)'./dynamicViscosity(nozzle.flow.T);
    
    % THE FOLLOWING 3 QUANTITIES NEED ESTIMATION
    nozzle.flow.hf = 0;
    nozzle.wallRecoveryFactor = 0;
    nozzle.flow.Cf = 0;
    
    % --- Extract wall data
    nozzle.Tw = nozzle.flow.T(Nj,:);
    % THE FOLLOWING QUANTITY NEEDS ESTIMATION, EXTERNAL WALL TEMP
    nozzle.Text = nozzle.Tw; % temporary fix

    % --- Record nozzle exit values
	nozzle.exit.M = avgDat(1,iMach);
	nozzle.exit.T = avgDat(1,iTem);
	nozzle.exit.U = (avgDat(1,iCons2)+avgDat(1,iCons3))/avgDat(1,iCons1);
	nozzle.exit.P = avgDat(1,iPres);
	nozzle.exit.Pstag = PstagExitAverage;
    nozzle.exit.Tstag = TstagExitAverage;    
    
    nozzle.PstagRatio = nozzle.exit.Pstag/nozzle.boundaryCdt.PtIn;
    nozzle.TstagRatio = nozzle.exit.Tstag/nozzle.boundaryCdt.TtIn;
    
    % --- Assign thrust & mass flow rate values
    nozzle.netThrust = thrust;
    nozzle.massFlowRate = mdot;
    
    % --- Calculate stresses
    % Stresses calculated assuming cylinder, nozzle length not constrained in 
    % thermal expansion
	%nozzle.hoopStress = prod([Pq(Nj,:) ;D(xq(1,:)')';1./(2*t(xq(1,:)')')]);
    nozzle.stress.hoop = Pq(Nj,:)'.*nozzle.geometry.D./(2*nozzle.wall.t);
    
    % Thermal stresses calculated assuming steady-state, give max tensile stress
    ri = nozzle.geometry.D/2; % inner radius
    ro = nozzle.geometry.D/2 + nozzle.wall.t; % outer radius
    nozzle.stress.thermal.radial = nozzle.wall.E*nozzle.wall.coeffThermalExpansion*(nozzle.Tw-nozzle.Text)/(2*(1-nozzle.wall.poissonRatio)).*(1./log(ro./ri)).*(1 - 2*ri.^2./(ro.^2 - ri.^2).*log(ro./ri));
    nozzle.stress.thermal.tangential = nozzle.stress.thermal.radial;

    nozzle.stress.maxPrincipal = nozzle.stress.hoop + nozzle.stress.thermal.tangential;
    nozzle.stress.principal = [nozzle.stress.maxPrincipal, nozzle.stress.thermal.radial, zeros(length(xPosition),1)]; 
    
    % --- Calculate cycles to failure Nf
    nozzle.Nf = estimateNf(nozzle.Tw,nozzle.stress.maxPrincipal,1);
	
    % --- Calc nozzle material volume
    % Volume calculation only works for spline parameterized nozzle geometry
    % and piecewise-linear parameterized nozzle wall thickness
    if(exist('pp','var')) % Exact volume for cubic spline parameterization
        nozzle.geometry.volume = wallVolume(pp,nozzle.wall);
    else % Approximate volume using trapezoidal integration
        xVolume = linspace(0,nozzle.geometry.length,500)';
        volumeIntegrand = pi*D(xVolume).*t(xVolume) + pi*t(xVolume).^2;
        nozzle.geometry.volume = (xVolume(2)-xVolume(1))*trapz(volumeIntegrand);
    end
      
	fprintf('\n -- Info: CFD results (averaged values at nozzle exit)\n');
	fprintf('           Mach   = %f\n', nozzle.exit.M );
	fprintf('           Temp   = %f K\n', nozzle.exit.T );
	fprintf('           Vel    = %f m/s\n', nozzle.exit.U );
	fprintf('           Pres   = %f Pa\n', nozzle.exit.P );
	fprintf('           Thrust = %f N\n', nozzle.netThrust );
	
end