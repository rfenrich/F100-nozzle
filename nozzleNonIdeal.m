function [ nozzle ] = nozzleNonIdeal( fluid, freestream, nozzle, err )
% Solve for flow along length of non-ideal nozzle given geometry, inlet
% stagnation temperature and pressure, and freestream temperature and
% pressure. Iterate for Cf and stagnation temperature. An ODE for M^2 is 
% solved given A, Cf, and Tstag. Pstag is found from mass conservation. T 
% and P are found from def'n of stag. temp. Density rho is found from ideal 
% gas law. 
%
% Returns M, density, pressure P, temperature T, stagnation 
% temp. Tstag, stagnation pressure Pstag, velocity U, Re, internal heat
% transfer coefficient hf, friction coefficient Cf, interior wall temp. Tw,
% exterior wall temp. Text, and approximate stress along length of nozzle.
%
% INPUTS:
% fluid = structure with fields: gam (ratio of specific heats) and R
%         (specific ideal gas constant)
% freestream = structure with fields: T and P
% nozzle = structure with fields: geometry.Ainlet2Athroat, 
%          geometry.Aexit2Athroat, geometry.length, geometry.shape,
%          geometry.xThroat,  geometry.xExit, inlet.Tstag, inlet.Pstag,
%          inlet.D, hInf
% err = structure defining error tolerances for iterations and solvers
%         included the following fields: err.betweenIterations.exitTemp,
%         err.solver.apparentThroatLocation, err.dMdxDenominator,
%         err.solver.M2relative, err.solver.M2absolute
%
% OUTPUTS:
% nozzle = modified input structure with additional fields including flow
% and specific geometry
%
% Rick Fenrich 7/17/15 modified 4/14/16

% ========================== GAS PROPERTIES ==============================
gam = fluid.gam;
R = fluid.R;

% Sutherland's law for dynamic viscosity:
dynamicViscosity = @(T) 1.716e-5*(T/273.15).^1.5*(273.15 + 110.4)./(T + 110.4); % kg/m*s

% ========================= NOZZLE PROPERTIES ============================
% Set inlet properties
inlet = nozzle.inlet;

% Calculate nozzle inlet, throat, and exit areas if they are not given:
if(~exist('nozzle.inlet.A','var'))
    nozzle.inlet.A = pi*inlet.D^2/4;
end
if(~exist('nozzle.throat.A','var'))
    nozzle.throat.A = nozzle.inlet.A/nozzle.geometry.Ainlet2Athroat;
end
if(~exist('nozzle.exit.A','var'))
    nozzle.exit.A = nozzle.geometry.Aexit2Athroat*nozzle.inlet.A/nozzle.geometry.Ainlet2Athroat;
end

nozzle.geometry.xApparentThroat = nozzle.geometry.xThroat; % initialize apparent throat location

% ========================== NOZZLE GEOMETRY =============================
[ A, dAdx, D, nozzle ] = nozzleParameterization( nozzle );
x = linspace(0,nozzle.geometry.length,200);
% ======================= NOZZLE WALL GEOMETRY ===========================
% Only piecewise-linear geometry is enabled thus far
t = @(x) piecewiseLinearGeometry(x,'t',nozzle.wall); % m, thickness of wall

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

% ==================== DETERMINE IDEAL NOZZLE FLOW =======================
pressureRatio = inlet.Pstag/freestream.P; % assumes no Pstag loss
nozzle.throat.Tstag = inlet.Tstag;
nozzle.throat.Pstag = inlet.Pstag;
nozzle.exit.Tstag = inlet.Tstag;
nozzle.exit.Pstag = inlet.Pstag;
[nozzle.status, nozzle.shock] = nozzleState( fluid, pressureRatio, nozzle.throat, nozzle.exit );
%fprintf('%s\n',nozzle.status);

% ======================= NOZZLE FRICTION & HEAT =========================

% Make a first guess at stagnation temperature and Cf along nozzle
dTstagdx = @(x) -6; % originally -6*x as first guess
Tstag = @(x) inlet.Tstag;
Pstag = @(x) inlet.Pstag;
Cf = @(x) 0.004;

xPositionOld = [0; nozzle.geometry.length];
flow.Tstag = [inlet.Tstag; inlet.Tstag];

converged = false;
maxIterations = 10; % max number of iterations to solve for Cf and Tstag
counter = 0; % used to count number of iterations
tolerance = err.betweenIterations.exitTemp; % tolerance for percent error in exit static temperature between iterations
Texit_old = 0; % saves previous exit static temperature

while ~converged
    
    counter = counter + 1;
    
    % ODE solver options
    options.RelTol = err.solver.M2relative;
    options.AbsTol = err.solver.M2absolute;  
    options.Events = @eventsFcn;
    
    % Find where M = 1
    nozzle.geometry.xApparentThroat = findApparentThroat(gam, nozzle, err, dAdx, A, Cf, D, dTstagdx, Tstag);
        
    if( strcmp(nozzle.status,'no flow') )
        
        error('Flow reversed in nozzle\n');
    
    % ------------------- INTEGRATE SUBSONIC FLOW ------------------------
    elseif( strcmp(nozzle.status,'subsonic') ) % subsonic flow
            
        if(inlet.Mach)
            
            fprintf('using given inlet Mach number\n');
            M0 = inlet.Mach;
            
            dM2dxPost = @(x, M2) (2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(x)./A(x) + 2*gam*M2*Cf(x)./D(x) + (1+gam*M2)*dTstagdx(x)./(2*Tstag(x)));
            
            % Solve using 4th-order Runge-Kutta method
            [xPosition,M2,TE,YE,IE] = ode45(dM2dxPost,[0 nozzle.geometry.length],M0.^2,options);
            
            if(IE == 1) % terminated when reached M = 0.9999
                error('given inlet Mach number too high');
            else
                checkIntegration(TE,YE,IE);
            end

        else
            
            %M0 = findMach(fluid, nozzle.geometry.xApparentThroat, err, nozzle.inlet, Pstag, A, Tstag);
 
            if(nozzle.geometry.xApparentThroat == nozzle.geometry.length) % only need to integrate forward

                dM2dxPrior = @(x, M2) -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(nozzle.geometry.xApparentThroat-x)./A(nozzle.geometry.xApparentThroat-x) + 2*gam*M2*Cf(nozzle.geometry.xApparentThroat-x)./D(nozzle.geometry.xApparentThroat-x) + (1+gam*M2)*dTstagdx(nozzle.geometry.xApparentThroat-x)./(2*Tstag(nozzle.geometry.xApparentThroat-x)));

                % Set starting conditions for integrations
                LowerM = 0.9999; % start integration at this Mach number

                % Solve using 4th-order Runge-Kutta method
                [xPositionPrior,M2Prior,TE,YE,IE] = ode45(dM2dxPrior,[0 nozzle.geometry.xApparentThroat],LowerM.^2,options);
                checkIntegration(TE,YE,IE);
                
                % Combine both parts of problem
                M2 = flipud(M2Prior); % contains Mach^2
                xPosition = -flipud(xPositionPrior)+nozzle.geometry.xApparentThroat;

            else % integrate forwards and backwards from throat

                % Split problem into before and after nozzle throat, solve for d(M^2)/dx
                dM2dxPrior = @(x, M2) -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(nozzle.geometry.xApparentThroat-x)./A(nozzle.geometry.xApparentThroat-x) + 2*gam*M2*Cf(nozzle.geometry.xApparentThroat-x)./D(nozzle.geometry.xApparentThroat-x) + (1+gam*M2)*dTstagdx(nozzle.geometry.xApparentThroat-x)./(2*Tstag(nozzle.geometry.xApparentThroat-x)));
                dM2dxPost = @(x, M2) (2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(x+nozzle.geometry.xApparentThroat)./A(x+nozzle.geometry.xApparentThroat) + 2*gam*M2*Cf(x+nozzle.geometry.xApparentThroat)./D(x+nozzle.geometry.xApparentThroat) + (1+gam*M2)*dTstagdx(x+nozzle.geometry.xApparentThroat)./(2*Tstag(x+nozzle.geometry.xApparentThroat)));

                % Set starting conditions for integrations
                UpperM = 0.9999; % start integration at this Mach number for aft portion of nozzle
                LowerM = 0.9999; % start integration at this Mach number for fore portion of nozzle
                dx = 1e-5; % 1e-5 for 0.9999 to 1.0001 or 1e-4 for 0.999 to 1.001       

                % Solve using 4th-order Runge-Kutta method
                [xPositionPost,M2Post,TE,YE,IE] = ode45(dM2dxPost,[dx/2 nozzle.geometry.length - nozzle.geometry.xApparentThroat],UpperM.^2,options);
                checkIntegration(TE,YE,IE);
                [xPositionPrior,M2Prior,TE,YE,IE] = ode45(dM2dxPrior,[dx/2 nozzle.geometry.xApparentThroat],LowerM.^2,options);
                checkIntegration(TE,YE,IE);
                
                % Combine both parts of problem
                M2 = [flipud(M2Prior); M2Post]; % contains Mach^2
                xPosition = [-flipud(xPositionPrior)+nozzle.geometry.xApparentThroat; xPositionPost + nozzle.geometry.xApparentThroat];

            end
        
        end                
            
    % ------------------- INTEGRATE OVER SHOCK ------------------------        
    elseif( strcmp(nozzle.status,'shock') ) % shock present 
            
        if(nozzle.geometry.xApparentThroat == nozzle.geometry.length) % only need to integrate forward
            
            dM2dxPrior = @(x, M2) -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(nozzle.geometry.xApparentThroat-x)./A(nozzle.geometry.xApparentThroat-x) + 2*gam*M2*Cf(nozzle.geometry.xApparentThroat-x)./D(nozzle.geometry.xApparentThroat-x) + (1+gam*M2)*dTstagdx(nozzle.geometry.xApparentThroat-x)./(2*Tstag(nozzle.geometry.xApparentThroat-x)));
            
            % Set starting conditions for integrations
            LowerM = 0.9999; % start integration at this Mach number
 
            % Solve using 4th-order Runge-Kutta method
            [xPositionPrior,M2Prior,TE,YE,IE] = ode45(dM2dxPrior,[0 nozzle.geometry.xApparentThroat],LowerM.^2,options);
            checkIntegration(TE,YE,IE);
            
            % Combine both parts of problem
            M2 = flipud(M2Prior); % contains Mach^2
            xPosition = -flipud(xPositionPrior)+nozzle.geometry.xApparentThroat;
            
        else % integrate forwards and backwards from throat

            % Split problem into before and after nozzle throat, solve for d(M^2)/dx
            dM2dxPost = @(x, M2) (2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(x+nozzle.geometry.xApparentThroat)./A(x+nozzle.geometry.xApparentThroat) + 2*gam*M2*Cf(x+nozzle.geometry.xApparentThroat)./D(x+nozzle.geometry.xApparentThroat) + (1+gam*M2)*dTstagdx(x+nozzle.geometry.xApparentThroat)./(2*Tstag(x+nozzle.geometry.xApparentThroat)));
            dM2dxPrior = @(x, M2) -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(nozzle.geometry.xApparentThroat-x)./A(nozzle.geometry.xApparentThroat-x) + 2*gam*M2*Cf(nozzle.geometry.xApparentThroat-x)./D(nozzle.geometry.xApparentThroat-x) + (1+gam*M2)*dTstagdx(nozzle.geometry.xApparentThroat-x)./(2*Tstag(nozzle.geometry.xApparentThroat-x)));

            % Set starting conditions for integrations
            UpperM = 1.0001; % start integration at this Mach number for aft portion of nozzle
            LowerM = 0.9999; % start integration at this Mach number for fore portion of nozzle
            dx = 1e-5; % 1e-5 for 0.9999 to 1.0001 or 1e-4 for 0.999 to 1.001       

            % Solve using 4th-order Runge-Kutta method
            [xPositionPrior,M2Prior,TE,YE,IE] = ode45(dM2dxPrior,[dx/2 nozzle.geometry.xApparentThroat],LowerM.^2,options);
            checkIntegration(TE,YE,IE);
            [xPositionPost,M2Post,TE,YE,IE] = ode45(dM2dxPost,[dx/2 nozzle.geometry.length - nozzle.geometry.xApparentThroat],UpperM.^2,options);
            checkIntegration(TE,YE,IE);
            
            % Generate shock and restart integration
            if(IE == 4) % Shock found
                % Calculate Mach number after shock
                xShockStart = TE;
                nozzle.shock.x = nozzle.geometry.xApparentThroat + TE;
                nozzle.shock.index = length(xPositionPrior) + length(xPositionPost) + 1;
                [xPositionPostShock,M2PostShock,TE,YE,IE] = ode45(dM2dxPost,[xShockStart nozzle.geometry.length - nozzle.geometry.xApparentThroat],nozzle.shock.Mpost^2,options);      
                checkIntegration(TE,YE,IE);
            else
                error('shock not found');
            end
            
            % Combine both parts of problem
            M2 = [flipud(M2Prior); M2Post; M2PostShock]; % contains Mach^2
            xPosition = [-flipud(xPositionPrior)+nozzle.geometry.xApparentThroat; xPositionPost + nozzle.geometry.xApparentThroat; xPositionPostShock + nozzle.geometry.xApparentThroat];
        
        end

    % ------------------- INTEGRATE SUPERSONIC FLOW ----------------------
    else % supersonic flow
        
        if(nozzle.geometry.xApparentThroat == nozzle.geometry.length) % only need to integrate forward
            
            dM2dxPrior = @(x, M2) -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(nozzle.geometry.xApparentThroat-x)./A(nozzle.geometry.xApparentThroat-x) + 2*gam*M2*Cf(nozzle.geometry.xApparentThroat-x)./D(nozzle.geometry.xApparentThroat-x) + (1+gam*M2)*dTstagdx(nozzle.geometry.xApparentThroat-x)./(2*Tstag(nozzle.geometry.xApparentThroat-x)));
            
            % Set starting conditions for integrations
            LowerM = 0.9999; % start integration at this Mach number
 
            % Solve using 4th-order Runge-Kutta method
            [xPositionPrior,M2Prior,TE,YE,IE] = ode45(dM2dxPrior,[0 nozzle.geometry.xApparentThroat],LowerM.^2,options);
            checkIntegration(TE,YE,IE);
            
            % Combine both parts of problem
            M2 = flipud(M2Prior); % contains Mach^2
            xPosition = -flipud(xPositionPrior)+nozzle.geometry.xApparentThroat;
            
        else % integrate forwards and backwards from throat

            % Split problem into before and after nozzle throat, solve for d(M^2)/dx
            dM2dxPrior = @(x, M2) -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(nozzle.geometry.xApparentThroat-x)./A(nozzle.geometry.xApparentThroat-x) + 2*gam*M2*Cf(nozzle.geometry.xApparentThroat-x)./D(nozzle.geometry.xApparentThroat-x) + (1+gam*M2)*dTstagdx(nozzle.geometry.xApparentThroat-x)./(2*Tstag(nozzle.geometry.xApparentThroat-x)));
            dM2dxPost = @(x, M2) (2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(x+nozzle.geometry.xApparentThroat)./A(x+nozzle.geometry.xApparentThroat) + 2*gam*M2*Cf(x+nozzle.geometry.xApparentThroat)./D(x+nozzle.geometry.xApparentThroat) + (1+gam*M2)*dTstagdx(x+nozzle.geometry.xApparentThroat)./(2*Tstag(x+nozzle.geometry.xApparentThroat)));

            % Set starting conditions for integrations
            UpperM = 1.0001; % start integration at this Mach number for aft portion of nozzle
            LowerM = 0.9999; % start integration at this Mach number for fore portion of nozzle
            dx = 1e-4; % 1e-5 for 0.9999 to 1.0001 or 1e-4 for 0.999 to 1.001

            % Determine dMdx across M = 1
    %         xBefore = nozzle.geometry.xApparentThroat - dx/2;
    %         xAfter =  nozzle.geometry.xApparentThroat + dx/2;
    %         dMdxBefore = (1 + (gam-1)/2)*LowerM^3/(1-LowerM^2)*(-dAdx(xBefore)/A(xBefore) + 2*gam*LowerM^2*Cf(xBefore)/D(xBefore) + (1+gam*LowerM^2)*dTstagdx(xBefore)/(2*Tstag(xBefore)));
    %         dMdxAfter = (1 + (gam-1)/2)*UpperM^3/(1-UpperM^2)*(-dAdx(xAfter)/A(xAfter) + 2*gam*UpperM^2*Cf(xAfter)/D(xAfter) + (1+gam*UpperM^2)*dTstagdx(xAfter)/(2*Tstag(xAfter)));
    %         dMdx = (abs(dMdxBefore)+abs(dMdxAfter))/2;
    %         dx2 = (UpperM-LowerM)/dMdx;
    %         fprintf('dx: %f\n',dx2);
    %         dx = dx2;        

            % Solve using 4th-order Runge-Kutta method
            [xPositionPost,M2Post,TE,YE,IE] = ode45(dM2dxPost,[dx/2 nozzle.geometry.length - nozzle.geometry.xApparentThroat],UpperM.^2,options);
            checkIntegration(TE,YE,IE);
            [xPositionPrior,M2Prior,TE,YE,IE] = ode45(dM2dxPrior,[dx/2 nozzle.geometry.xApparentThroat],LowerM.^2,options);
            checkIntegration(TE,YE,IE);
            
            % Combine both parts of problem
            M2 = [flipud(M2Prior); M2Post]; % contains Mach^2
            xPosition = [-flipud(xPositionPrior)+nozzle.geometry.xApparentThroat; xPositionPost + nozzle.geometry.xApparentThroat];
        
        end
        
    end
    
    % Throw exception if M^2 is negative for whatever reason
    if (sum(M2<0) > 0) % i.e. if any M2 < 0, implying Mach number is imaginary
        error('! Imaginary Mach number calculated')
    end

    % Plot for testing
%     hold on
%     plot(xPosition,sqrt(M2));
%     drawnow;
    
    % ======================= CALC OTHER PROPERTIES ==========================

    flow.M = sqrt(M2); % Mach number
    flow.Tstag = interpLinear(xPositionOld,flow.Tstag,xPosition);
    %flow.Tstag = interp1(xPositionOld,flow.Tstag,xPosition,'linear');
    flow.T = flow.Tstag./(1 + (gam-1)*M2/2); % static temperature from stag. temp. definition
    flow.Pstag = inlet.Pstag*(A(0)./A(xPosition)).*(AreaMachFunc(gam,flow.M(1))./AreaMachFunc(gam,flow.M)).*sqrt(flow.Tstag/flow.Tstag(1)); % stagnation pressure from mass conservation        
    flow.P = flow.Pstag./(1 + (gam-1)*M2/2).^(gam/(gam-1)); % static pressure from stag. press. definition
    flow.density = flow.P./(R*flow.T); % density
    flow.U = flow.M.*sqrt(gam*R*flow.T); % velocity
    flow.Re = flow.density.*flow.U.*D(xPosition)./dynamicViscosity(flow.T); % Reynolds number from definition
		
    % =================== RECALCULATE FRICTION & HEAT ========================

    % Heat transfer
    T = @(x) interpLinear(xPosition,flow.T,x);
    %T = @(x) interp1(xPosition,flow.T,x,'linear');
    flow.hf = Pr(flow.T).^(2/3).*flow.density.*Cp(flow.T).*flow.U.*Cf(xPosition)/2; % heat transfer coefficient to interior nozzle wall, estimated using Chilton-Colburn analogy
    
    % Redefine stagnation temperature distribution
    TstagXIntegrand = 4./(Cp(flow.T).*flow.density.*flow.U.*D(xPosition).*(1./flow.hf + t(xPosition)/nozzle.wall.k + 1/nozzle.hInf));
    TstagXIntegral = cumtrapz(xPosition,TstagXIntegrand);
    flow.Tstag = freestream.T*(1 - exp(-TstagXIntegral)) + inlet.Tstag*exp(-TstagXIntegral);
    Tstag = @(x) interpLinear(xPosition,flow.Tstag,x);
    Pstag = @(x) interpLinear(xPosition,flow.Pstag,x);
    %Tstag = @(x) interp1(xPosition,flow.Tstag,x,'linear');
    dTstagdxVal = (freestream.T - flow.Tstag)*4./(Cp(flow.T).*flow.density.*flow.U.*D(xPosition).*(1./flow.hf + t(xPosition)/nozzle.wall.k + 1/nozzle.hInf));
    dTstagdx = @(x) interpLinear(xPosition,dTstagdxVal,x);
    %dTstagdx = @(x) interp1(xPosition,dTstagdxVal,x,'linear');
    
    % Estimate interior wall temperature
    Qw = Cp(flow.T).*flow.density.*flow.U.*D(xPosition).*dTstagdxVal/4;
    nozzle.wall.Tinside = flow.Tstag + Qw./flow.hf; % wall temperature
    nozzle.wall.recoveryFactor = (nozzle.wall.Tinside./flow.T - 1)./((gam-1)*flow.M.^2/2);
    
    % Estimate exterior wall temperature
    nozzle.wall.Toutside = ones(length(xPosition),1)*freestream.T - Qw./nozzle.hInf; % nozzle exterior wall temp.
    
    % Redefine friction coefficient distribution (Sommer & Short's method)
    TPrimeRatio = 1 + 0.035*flow.M.^2 + 0.45*(nozzle.wall.Tinside./flow.T -1);
    RePrimeRatio = 1./(TPrimeRatio.*(TPrimeRatio).^1.5.*(1 + 110.4./flow.T)./(TPrimeRatio + 110.4./flow.T));
    CfIncomp = 0.074./flow.Re.^0.2;
    flow.Cf = CfIncomp./TPrimeRatio./RePrimeRatio.^0.2;
    Cf = @(x) interpLinear(xPosition,flow.Cf,x);
    %Cf = @(x) interp1(xPosition,flow.Cf,x,'linear');
    
    % Save old solution x position for next iteration
    xPositionOld = xPosition;
    
    % Recalculate nozzle status
    nozzle.throat.Tstag = Tstag(nozzle.geometry.xApparentThroat);
    nozzle.throat.Pstag = Pstag(nozzle.geometry.xApparentThroat);
    nozzle.exit.Tstag = flow.Tstag(end);
    if(nozzle.shock.M) % shock in nozzle
        % approximate calculation
        pressureRatio = flow.Pstag(nozzle.shock.index-1)/freestream.P;
        nozzle.exit.Pstag = flow.Pstag(nozzle.shock.index-1); % temp. fix
        [nozzle.status, nozzle.shock] = nozzleState( fluid, pressureRatio, nozzle.throat, nozzle.exit );
        nozzle.exit.Pstag = flow.Pstag(end);
    else % no shock in nozzle
        pressureRatio = flow.Pstag(end)/freestream.P; % assumes no Pstag loss
        nozzle.exit.Pstag = flow.Pstag(end);
        [nozzle.status, nozzle.shock] = nozzleState( fluid, pressureRatio, nozzle.throat, nozzle.exit );
    end
    %fprintf('%s\n',nozzle.status);
    
    fprintf('%i ',counter);
    
    if counter >= maxIterations
        fprintf('! Max iterations reached for non-ideal nozzle heat transfer & friction.\n');
        break;
    end
    
    % Check tolerance on static temperature at nozzle exit
    percentError = abs(T(nozzle.geometry.length)-Texit_old)/T(nozzle.geometry.length);
    %fprintf('%% Error: %e\n',percentError*100);
    Texit_old = T(nozzle.geometry.length);
    if(percentError < tolerance)
        converged = true;
        fprintf('iter to converge non-ideal nozzle heat xfer & friction calcs\n');
    end
    
end

% =================== ASSIGN FLOW DATA TO NOZZLE =========================

nozzle.flow = flow;
nozzle.PstagRatio = flow.Pstag(end)/flow.Pstag(1);
nozzle.TstagRatio = flow.Tstag(end)/flow.Tstag(1);

% Record nozzle exit values
nozzle.exit.M = nozzle.flow.M(end);
nozzle.exit.Pstag = flow.Pstag(end);
nozzle.exit.Tstag = flow.Tstag(end);
nozzle.exit.P = nozzle.flow.P(end);
nozzle.exit.T = nozzle.flow.T(end);
nozzle.exit.U = nozzle.flow.U(end);

% ========================== CALC GEOMETRY ===============================

nozzle.xPosition = xPosition;
nozzle.geometry.A = A(xPosition);
nozzle.geometry.dAdx = dAdx(xPosition);
nozzle.geometry.D = D(xPosition);
nozzle.wall.t = t(xPosition);
nozzle.geometry.maxSlope = max(nozzle.geometry.dAdx./pi./nozzle.geometry.D);
nozzle.geometry.minSlope = min(nozzle.geometry.dAdx./pi./nozzle.geometry.D);

% ==================== CALCULATE APPROX. THRUST ==========================
% Thrust is only approximate since engine is not taken into account.

% Estimate losses due to divergence (formula below is for conical nozzle)
nozzle.exit.angle = atan(dAdx(nozzle.geometry.length)/pi/D(nozzle.geometry.length))*180/pi;
nozzle.divergenceFactor = (1 + cosd(nozzle.exit.angle))/2;

nozzle.massFlowRate = massFlowRate(fluid,nozzle.inlet.Pstag,nozzle.inlet.A,nozzle.inlet.Tstag,nozzle.flow.M(1));
nozzle.netThrust = nozzle.divergenceFactor*nozzle.massFlowRate*(nozzle.exit.U - freestream.U) + (nozzle.exit.P - freestream.P)*nozzle.exit.A;
nozzle.grossThrust = nozzle.divergenceFactor*nozzle.massFlowRate*(nozzle.exit.U) + (nozzle.exit.P - freestream.P)*nozzle.exit.A;

% ========================== CALC STRESSES ===============================
% Stresses calculated assuming cylinder, nozzle length not constrained in 
% thermal expansion

nozzle.stress.hoop = flow.P.*D(xPosition)./(2*t(xPosition));

% Thermal stresses calculated assuming steady-state, give max tensile stress
ri = nozzle.geometry.D/2; % inner radius
ro = nozzle.geometry.D/2 + nozzle.wall.t; % outer radius
nozzle.stress.thermal.radial = nozzle.wall.E*nozzle.wall.coeffThermalExpansion*(nozzle.wall.Tinside-nozzle.wall.Toutside)/(2*(1-nozzle.wall.poissonRatio)).*(1./log(ro./ri)).*(1 - 2*ri.^2./(ro.^2 - ri.^2).*log(ro./ri));
nozzle.stress.thermal.tangential = nozzle.stress.thermal.radial;

% Estimate vonMises, even though not really valid for composites
%nozzle.stress.vonMises = sqrt( (nozzle.stress.hoop+nozzle.stress.thermal.tangential).^2 - nozzle.stress.hoop.*nozzle.stress.thermal.radial + nozzle.stress.thermal.radial.^2 );
nozzle.stress.maxPrincipal = nozzle.stress.hoop + nozzle.stress.thermal.tangential;
nozzle.stress.principal = [nozzle.stress.maxPrincipal, nozzle.stress.thermal.radial, zeros(length(xPosition),1)];

% ==================== CALC CYCLES TO FAILURE NF =========================

nozzle.Nf = estimateNf(nozzle.wall.Tinside,nozzle.stress.maxPrincipal,1);

% =================== CALC NOZZLE MATERIAL VOLUME ========================

% Volume calculation only works for spline parameterized nozzle geometry
% and piecewise-linear parameterized nozzle wall thickness
if(exist('pp','var')) % Exact volume for cubic spline parameterization
    nozzle.geometry.volume = wallVolume(pp,nozzle.wall);
else % Approximate volume using trapezoidal integration
    nozzle.geometry.volume = wallVolume(nozzle.geometry.length,D,t,'integrate');
end

% ========================== PLOT GEOMETRY ===============================
%plot(nozzle.xPosition,nozzle.geometry.D/2)
%axis equal
%drawnow

    function [ value, isTerminal, direction ] = eventsFcn(x, M2)

        % Event function 1: Mach becomes larger than 0.9999
        value(1) = sqrt(M2) - 0.9999;
        isTerminal(1) = 0; 
        direction(1) = 1; % only flag an event if M2 is increasing

        %Event function 2: Mach becomes smaller than 1.0001
        value(2) = sqrt(M2) - 1.0001;
        isTerminal(2) = 0;
        direction(2) = -1; % only flag an event if M2 is decreasing

        % Event function 3: Mach is 1
        value(3) = sqrt(M2) - 1.00;
        isTerminal(3) = 1;
        direction(3) = 0; % flag if any

        % Event function 4: Mach shock number reached
        value(4) = sqrt(M2) - nozzle.shock.M;
        isTerminal(4) = 1;
        direction(4) = 1; % only flag an event if M2 is increasing

    end

end

function [result] = AreaMachFunc(g,M)
    result = ((g+1)/2)^((g+1)/(2*(g-1)))*M./(1+(g-1)*M.^2/2).^((g+1)/(2*(g-1)));
end

function [result] = massFlowRate(fluid,Pt,Area,Tt,M)
    gam = fluid.gam;
    R = fluid.R;
    result = (gam/((gam+1)/2)^((gam+1)/(2*(gam-1))))*Pt*Area*AreaMachFunc(gam,M)/sqrt(gam*R*Tt);
end

function [xApparentThroat] = findApparentThroat(gam, nozzle, err, dAdx, A, Cf, D, dTstagdx, Tstag)

    % Solve for location where M = 1 (location of apparent throat)
    options = optimset('TolFun',err.solver.apparentThroatLocation);
    
    dMdxCoeffFunc = @(x) -dAdx(x)./A(x) + 2*gam*Cf(x)./D(x) + (1+gam)*dTstagdx(x)./(2*Tstag(x));
    
    % Find sign changes in dMdxCoeffFunc
    xFind = linspace(0,nozzle.geometry.length-1e-4,200)';
    coeffFind = dMdxCoeffFunc(xFind);
    coeffFind(coeffFind>0) = 1;
    coeffFind(coeffFind<=0) = 0;
    signChangeLocations = find(diff(coeffFind)~=0);
    if(isempty(signChangeLocations))
        throatGuess = nozzle.geometry.xThroat;
    else
        [minVal,minInd] = min(A(xFind(signChangeLocations)));
        throatGuess = xFind(signChangeLocations(minInd));
    
        % Check to make sure each following possible throat is far enough away
        ind = signChangeLocations(minInd);
        for ii = minInd+1:length(signChangeLocations)

            currentInd = signChangeLocations(ii);
            dx = xFind(currentInd) - xFind(ind);

            dAdxbar = dAdx(xFind(ind));
            dAdxbarEst = (A(xFind(currentInd)) - A(xFind(ind)))/dx;

            dTstagdxbar = dTstagdx(xFind(currentInd));
            Abar = A(xFind(ind));
            Tstagbar = Tstag(xFind(ind));
            Cfbar = Cf(xFind(ind));
            Dbar = D(xFind(ind));

            RHS = dTstagdxbar*Abar*(1+gam)/(2*Tstagbar) + 2*gam*Abar*Cfbar/Dbar;

            if( dAdxbarEst <= RHS )
                throatGuess = xFind(currentInd);
            end

        end
    
    end
    
    %xApparentThroat = fzero(dMdxCoeffFunc,nozzle.geometry.xThroat,options);
    xApparentThroat = fzero(dMdxCoeffFunc,throatGuess,options);
    if( isnan(xApparentThroat) || xApparentThroat > nozzle.geometry.length)
        xApparentThroat = nozzle.geometry.length;
        fprintf('Mach = 1 at exit\n');
    end
    
    if(D(xApparentThroat) > D(nozzle.geometry.length))
        xApparentThroat = nozzle.geometry.length;
    end
    
end

function [inletMach] = findMach(fluid, xApparentThroat, err, nozzleConditions, Pstag, A, Tstag)
    % Equates mass flow rate at apparent throat and mass flow rate at
    % desired location which has known conditions: Pstag, Tstag, and area A
    
    gam = fluid.gam;
    R = fluid.R;
    
    options = optimset('TolFun',err.solver.apparentThroatLocation);

    % Solve for mdot
    mdot = massFlowRate(fluid,Pstag(xApparentThroat),A(xApparentThroat),Tstag(xApparentThroat),1);
    
    % Solve for inlet Mach
    inletMachFcn = @(Mach) mdot - massFlowRate(fluid,nozzleConditions.Pstag,nozzleConditions.A,nozzleConditions.Tstag,Mach);
    inletMach = fzero(inletMachFcn,0.4,options);
    fprintf('inletMach: %f\n',inletMach);
    
end

function [] = checkIntegration(TE,YE,IE)

    if(~isempty(IE))
        
        for ii = 1:length(IE)
            
            if(IE(ii) == 1)

                fprintf('\n\n !!! Mach exceeds 0.9999 at x = %f',TE(ii));
                %error('Mach exceeds 0.9999 at x = %f',TE(ii));
                
            elseif(IE(ii) == 2)
                
                fprintf('\n\n !!! Mach falls below 1.0001 at x = %f',TE(ii));
                %error('Mach falls below 1.0001 at x = %f',TE(ii));
                
            elseif(IE(ii) == 3)
                
                fprintf('\n\n !!! Mach reaches 1 at x = %f',TE(ii));
                %error('Mach reaches 1 at x = %f',TE(ii));
                
            elseif(IE(ii) == 4)
                
                % no error, do nothing
                
            else
                
                error('Only 4 event functions are specified for ODE integration');
                
            end
            
        end
        
    end

end
