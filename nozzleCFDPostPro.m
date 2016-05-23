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

    x = dat(:,ix);
    y = dat(:,iy);
    v = dat(:,iTem);

    rhoU =  dat(:,iCons2);
    rho  =  dat(:,iCons1);
    P    =  dat(:,iPres);
    M    =  dat(:,iMach);
    T    =  dat(:,iTem);

    % --------------------------------------------------------
    % ---    Interpolate solution to structured grid
    % --------------------------------------------------------

    Ni=200;
    Nj=10;

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
    %
    % /!\ Not in use anymore, as it is not possible to maintain an inner marker aligned with the exit
    %     because of the length varies during the optimization
    %
    %idx = find( dat(:,ix) == xextract & dat(:,iy) < D(nozzle.geometry.length)/2 + 1e-6 );
    %DatLin = dat(idx,:);
    %DatLin = sortrows(DatLin,iy);
    %
    %NbvLin = size(idx,1);
    %
    %if ( NbvLin < 2 ) 
    %	error('  ## ERROR nozzleCFDPostPro : No extraction data was found');
    %end 
    %
    %avgDat=zeros(1,size(DatLin,2));
    %
    %yTotal = max(DatLin(:,iy));
    %for i=2:NbvLin
    %	dy = DatLin(i,iy)-DatLin(i-1,iy);
    %	avgDat(1,:) = avgDat(1,:)+dy*DatLin(i,:);
    %  
    %  rhoU = DatLin(i,iCons2);
    %	rho  = DatLin(i,iCons1);
    %	U    = DatLin(i,iCons2)/DatLin(i,iCons1);
    %	U0   = freestream.U;
    %	P    = DatLin(i,iPres);
    %	P0   = freestream.P;
    %  M    = DatLin(i,iMach);
    %  T    = DatLin(i,iTem);
    %     
    %  % --- Compute stagnation properties at exit
    %  PstagExitAverage = PstagExitAverage + (dy/yTotal)*P*(1 + (fluid.gam-1)*M^2/2)^(fluid.gam/(fluid.gam-1));
    %  TstagExitAverage = TstagExitAverage + (dy/yTotal)*T*(1 + (fluid.gam-1)*M^2/2); 
    %	
    %	% --- Compute thrust
    %	% T = 2PI * Int_{0}^{R} (rho U ( U - U0) + P - Po ) r dr
    %	thrust = thrust + dy*(rhoU*(U-U0)+P-P0);
    %  mdot = mdot + dy*rhoU;
    %end

    %hei = DatLin(NbvLin,iy)-DatLin(1,iy);
    %avgDat(1,:) = avgDat(1,:)/hei;

    % -------------------------------------------------------------
    % ---- Extract solution at exit without using the inner edges
    % -------------------------------------------------------------

    rhoV =  dat(:,iCons3);
    rhoU =  dat(:,iCons2);
    rho  =  dat(:,iCons1);
    P    =  dat(:,iPres);
    M    =  dat(:,iMach);
    T    =  dat(:,iTem);

    Nii = 3;	
    Njj = int16(0.5*D(nozzle.geometry.length)/nozzle.sizWal)+1;
    Njj = min(max(Njj,4),200);

    xqq = repmat(linspace(0,nozzle.geometry.length,Nii),Njj,1);
    yqq = zeros(Njj,Nii);

    for i=1:size(xqq,2)
        yqq(:,i) = linspace(0, 0.5*D(xqq(1,i)), Njj);
    end

    U0   = freestream.U;
    P0   = freestream.P;

    rhoUqq =  griddata(x,y,rhoU,xqq,yqq);
    rhoVqq =  griddata(x,y,rhoV,xqq,yqq);
    rhoqq  =  griddata(x,y,rho ,xqq,yqq);
    Pqq    =  griddata(x,y,P   ,xqq,yqq);
    Mqq    =  griddata(x,y,M   ,xqq,yqq);
    Tqq    =  griddata(x,y,T   ,xqq,yqq);

    yTotal = yqq(Njj,Nii);
    dy     = yTotal/double(Njj-1);

    avgDat=zeros(1,size(dat,2));

    for j=2:Njj

        avgDat(1,iCons3) = avgDat(1,iCons3)+dy*rhoVqq(j,Nii);
        avgDat(1,iCons2) = avgDat(1,iCons2)+dy*rhoUqq(j,Nii);
        avgDat(1,iCons1) = avgDat(1,iCons1)+dy*rhoqq(j,Nii);
        avgDat(1,iPres ) = avgDat(1,iPres )+dy*Pqq(j,Nii);
        avgDat(1,iMach ) = avgDat(1,iMach )+dy*Mqq(j,Nii);
        avgDat(1,iTem  ) = avgDat(1,iTem  )+dy*Tqq(j,Nii);

        rhoU = rhoUqq (j,Nii);
        rho  = rhoqq  (j,Nii);
        P    = Pqq    (j,Nii);
        M    = Mqq    (j,Nii);
        T    = Tqq    (j,Nii);
        U    = rhoU/rho;

        % --- Compute stagnation properties at exit
        PstagExitAverage = PstagExitAverage + (dy/yTotal)*P*(1 + (fluid.gam-1)*M^2/2)^(fluid.gam/(fluid.gam-1));
        TstagExitAverage = TstagExitAverage + (dy/yTotal)*T*(1 + (fluid.gam-1)*M^2/2); 

        % --- Compute thrust
        % T = 2PI * Int_{0}^{R} (rho U ( U - U0) + P - Po ) r dr

        thrust = thrust + dy*(rhoU*(U-U0)+P-P0);
        mdot = mdot + dy*rhoU;
        
    end

    avgDat(1,:) = avgDat(1,:)/yTotal;

    %WriteGMFStruc ('toto.mesh', xqq, yqq, Pqq );

    % --- Extract pressure, temp, density, and velocity onto a grid
    P = dat(:,iPres);
    Pq = griddata(x,y,P,xq,yq);
    T = dat(:,iTem);
    Tq = griddata(x,y,T,xq,yq);
    rho = dat(:,iCons1);
    rhoq = griddata(x,y,rho,xq,yq);
    U = dat(:,iCons2)./dat(:,iCons1);
    Uq = griddata(x,y,U,xq,yq);

    %WriteGMFStruc ('tata.mesh', xq, yq, Pq )

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
    nozzle.flow.Re = rhoq.*abs(Uq).*repmat(nozzle.geometry.D,1,Nj)'./dynamicViscosity(nozzle.flow.T);

    % ========================= THERMAL PROPERTIES ===========================
    % Uncomment the following if you want Pr, conductivity k, and Cp to change
    % with temperature:
    air.temp = [175 200 225 250 275 300 325 350 375 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1500];
    air.Pr = [0.744 0.736 0.728 0.72 0.713 0.707 0.701 0.697 0.692 0.688 0.684 0.68 0.68 0.68 0.682 0.684 0.687 0.69 0.693 0.696 0.699 0.702 0.704 0.707 0.709 0.711 0.713 0.715 0.717 0.719 0.722];
    air.k = 0.01*[1.593 1.809 2.02 2.227 2.428 2.624 2.816 3.003 3.186 3.365 3.71 4.041 4.357 4.661 4.954 5.236 5.509 5.774 6.03 6.276 6.52 6.754 6.985 7.209 7.427 7.64 7.849 8.054 8.253 8.45 8.831];
    air.Cp = [1002.3 1002.5 1002.7 1003.1 1003.8 1004.9 1006.3 1008.2 1010.6 1013.5 1020.6 1029.5 1039.8 1051.1 1062.9 1075.0 1087.0 1098.7 1110.1 1120.9 1131.3 1141.1 1150.2 1158.9 1167.0 1174.6 1181.7 1188.4 1194.6 1200.5 1211.2];
    Pr = @(T) interpLinear(air.temp,air.Pr,T); % Prandtl number of air
    kf = @(T) interpLinear(air.temp,air.k,T); % thermal conductivity of air
    Cp = @(T) interpLinear(air.temp,air.Cp,T); % specific heat of air

    % Assume average values for Pr, thermal conductivity k, and Cp:
    %Pr = @(T) 0.7;
    %kf = @(T) 0.037;
    %Cp = @(T) 1035;

    % Estimate Cf
    TPrimeRatio = 1 + 0.035*nozzle.flow.M(Nj,:).^2 + 0.45*(nozzle.flow.Tstag(Nj,:)./nozzle.flow.T(Nj,:) -1);
    RePrimeRatio = 1./(TPrimeRatio.*(TPrimeRatio).^1.5.*(1 + 110.4./nozzle.flow.T(Nj,:))./(TPrimeRatio + 110.4./nozzle.flow.T(Nj,:)));
    CfIncomp = 0.074./nozzle.flow.Re(Nj,:).^0.2;
    nozzle.flow.Cf = CfIncomp./TPrimeRatio./RePrimeRatio.^0.2;
  
    % Estimate wall data (no iterations; this estimate is conservative)
    nozzle.flow.hf = 0.5*nozzle.flow.density(Nj,:)'.*abs(nozzle.flow.U(Nj,:)').*Cp(nozzle.flow.T(Nj,:)).*Pr(nozzle.flow.T(Nj,:)).^(2/3).*nozzle.flow.Cf';
    Rtotal = 1./nozzle.flow.hf + nozzle.wall.t./nozzle.wall.k + 1./nozzle.hInf; % total thermal resistance
    Qw = (freestream.T - nozzle.flow.Tstag(Nj,:)')./Rtotal; % heat passing through wall
    nozzle.wall.Tinside = nozzle.flow.Tstag(Nj,:)' + Qw./nozzle.flow.hf; % nozzle interior wall temp.
    nozzle.wall.Toutside = freestream.T - Qw./nozzle.hInf; % nozzle exterior wall temp.

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
    nozzle.stress.hoop = Pq(Nj,:)'.*nozzle.geometry.D./(2*nozzle.wall.t);

    % Thermal stresses calculated assuming steady-state, give max tensile stress
    ri = nozzle.geometry.D/2; % inner radius
    ro = nozzle.geometry.D/2 + nozzle.wall.t; % outer radius
    nozzle.stress.thermal.radial = nozzle.wall.E*nozzle.wall.coeffThermalExpansion*(nozzle.wall.Tinside-nozzle.wall.Toutside)/(2*(1-nozzle.wall.poissonRatio)).*(1./log(ro./ri)).*(1 - 2*ri.^2./(ro.^2 - ri.^2).*log(ro./ri));
    nozzle.stress.thermal.tangential = nozzle.stress.thermal.radial;

    nozzle.stress.maxPrincipal = nozzle.stress.hoop + nozzle.stress.thermal.tangential;
    nozzle.stress.principal = [nozzle.stress.maxPrincipal, nozzle.stress.thermal.radial, zeros(length(nozzle.xPosition),1)]; 

    % --- Calculate cycles to failure Nf
    nozzle.Nf = estimateNf(nozzle.wall.Tinside,nozzle.stress.maxPrincipal,1);
	  
	fprintf('\n -- Info: CFD results (averaged values at nozzle exit)\n');
	fprintf('           Mach   = %f\n', nozzle.exit.M );
	fprintf('           Temp   = %f K\n', nozzle.exit.T );
	fprintf('           Vel    = %f m/s\n', nozzle.exit.U );
	fprintf('           Pres   = %f Pa\n', nozzle.exit.P );
	fprintf('           Thrust = %f N\n', nozzle.netThrust );
	fprintf('           Pstag  = %f N\n', PstagExitAverage );
	fprintf('           Tstag  = %f N\n', TstagExitAverage );
	
end