function [ nozzle ] = nozzleNonIdeal( fluid, freestream, nozzle, error )
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
% exterior wall temp. Text, and hoopStress along length of nozzle.
%
% INPUTS:
% fluid = structure with fields: gam (ratio of specific heats) and R
%         (specific ideal gas constant)
% freestream = structure with fields: T and P
% nozzle = structure with fields: geometry.Ainlet2Athroat, 
%          geometry.Aexit2Athroat, geometry.length, geometry.shape,
%          geometry.xThroat,  geometry.xExit, inlet.Tstag, inlet.Pstag,
%          inlet.D, hInf
% error = structure defining error tolerances for iterations and solvers
%         included the following fields: error.betweenIterations.exitTemp,
%         error.solver.apparentThroatLocation, error.dMdxDenominator,
%         error.solver.M2relative, error.solver.M2absolute
%
% OUTPUTS:
% nozzle = modified input structure with additional fields including flow
% and specific geometry
%
% Rick Fenrich 7/17/15 modified 10/16/15

% ========================== GAS PROPERTIES ==============================
gam = fluid.gam;
R = fluid.R;

% Area-Mach function from mass 1-D mass conservation equations:
AreaMachFunc = @(g,M) ((g+1)/2)^((g+1)/(2*(g-1)))*M./(1+(g-1)*M.^2/2).^((g+1)/(2*(g-1)));
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

% Calculate pressure ratio which determines state of nozzle:
pressureRatio = inlet.Pstag/freestream.P;

% ========================== NOZZLE GEOMETRY =============================

if(strcmp(nozzle.geometry.shape,'spline'))
    % set up spline for nozzle
    if(ischar(nozzle.geometry.spline.seed)) % seed shape is given
        % Set control points for splines (xNode), value of function at control
        % point (yNode), and slopes at start and end of spline (slopes)
        xNode = nozzle.geometry.spline.breaks;
        if(max(xNode) > nozzle.geometry.length || min(xNode) < 0) % check user given location of control points
            error('Spline control point outside nozzle length domain');
        end
        yNode = nozzleGeometry(xNode,'D',nozzle.inlet.D,nozzle.geometry.length,nozzle.geometry.xThroat,nozzle.geometry.Ainlet2Athroat,nozzle.geometry.Aexit2Athroat,nozzle.geometry.spline.seed)/2;
        nozzle.geometry.spline.seed = [xNode, yNode];
    else
        % Extract control points for splines and their values from the
        % given array
        xNode = nozzle.geometry.spline.seed(:,1);
        yNode = nozzle.geometry.spline.seed(:,2);
        if(xNode(1) == 0 && xNode(end) == nozzle.geometry.length) % check user given control point values (yNode) match user given area ratios
            areaRatio = yNode(end)^2/yNode(1)^2;
            areaRatioTolerance = 1e-3;
            if(areaRatio > nozzle.geometry.Aexit2Athroat/nozzle.geometry.Ainlet2Athroat + areaRatioTolerance || areaRatio < nozzle.geometry.Aexit2Athroat/nozzle.geometry.Ainlet2Athroat - areaRatioTolerance)
               error('Spline control point values do not match given nozzle area ratios'); 
            end
        end
    end
    
    slopes = nozzle.geometry.spline.slopes;
    pp = spline(xNode,[slopes(1); yNode; slopes(2)]); % perform piecewise cubic spline interpolation
    
    % Adjust nozzle throat size/location information if it has changed
    [xThroat, yThroat] = splineGeometry(0, 'throat', pp);
    if(xThroat ~= nozzle.geometry.xThroat)
        fprintf('throat size/location changed with spline parameterization\n');
    end
    nozzle.geometry.xThroat = xThroat;
    nozzle.throat.A = pi*yThroat^2;
    nozzle.geometry.Ainlet2Athroat = nozzle.inlet.A/nozzle.throat.A;
    nozzle.geometry.Aexit2Athroat = nozzle.exit.A/nozzle.throat.A;

    % Make necessary functions for splined nozzle shape
    A = @(x) splineGeometry(x,'A',pp);
    dAdx = @(x) splineGeometry(x,'dAdx',pp);
    D = @(x) splineGeometry(x,'D',pp);

else % if nozzle shape is not a spline
    A = @(x) nozzleGeometry(x,'A',nozzle.inlet.D,nozzle.geometry.length,nozzle.geometry.xThroat,nozzle.geometry.Ainlet2Athroat,nozzle.geometry.Aexit2Athroat,nozzle.geometry.shape);
    dAdx = @(x) nozzleGeometry(x,'dAdx',nozzle.inlet.D,nozzle.geometry.length,nozzle.geometry.xThroat,nozzle.geometry.Ainlet2Athroat,nozzle.geometry.Aexit2Athroat,nozzle.geometry.shape);
    D = @(x) nozzleGeometry(x,'D',nozzle.inlet.D,nozzle.geometry.length,nozzle.geometry.xThroat,nozzle.geometry.Ainlet2Athroat,nozzle.geometry.Aexit2Athroat,nozzle.geometry.shape);
end

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
options.Display='none'; % used for fsolve
shockInNozzle = false;

exit.critical.Msubsonic = fsolve(@(x) AreaMachFunc(gam,x) - nozzle.throat.A/nozzle.exit.A, 0.5, options);
exit.critical.Msupersonic = fsolve(@(x) AreaMachFunc(gam,x) - nozzle.throat.A/nozzle.exit.A, 2, options);
exit.critical.PtRatioSubsonic = (1 + (gam-1)*exit.critical.Msubsonic^2/2)^(gam/(gam-1));
exit.critical.PtRatioSupersonic = (1 + (gam-1)*exit.critical.Msupersonic^2/2)^(gam/(gam-1));

MbehindShock = sqrt((1 + (gam-1)*exit.critical.Msupersonic^2/2)/(gam*exit.critical.Msupersonic^2 - (gam-1)/2));
exit.normalShock.PtRatio = (nozzle.exit.A/nozzle.throat.A)*((gam+1)/2)^((gam+1)/(2*(gam-1)))*MbehindShock*sqrt(1 + (gam-1)*MbehindShock^2/2);

deltaPtRatio = 0.05; % a tolerance for deciding whether fully expanded flow occurs

if (pressureRatio <= 1)
    nozzle.status = 'no flow';
    error('! Nozzle inlet pressure less than freestream pressure.');
elseif (pressureRatio < exit.critical.PtRatioSubsonic)
    %fprintf('ideal: Subsonic flow throughout\n');
    nozzle.status = 'ideal: subsonic';
    exit.M = sqrt(2/(gam-1))*sqrt(0.999*pressureRatio^((gam-1)/gam) - 1); % first guess at exit.M
elseif (pressureRatio < exit.normalShock.PtRatio)
    %fprintf('ideal: Shock in nozzle\n');
    nozzle.status = 'ideal: shock in nozzle';
    exit.M = fsolve( @(x) ((gam+1)/2)^((gam+1)/(2*(gam-1)))*x*sqrt(1 + (gam-1)*x^2/2) - pressureRatio*nozzle.throat.A/nozzle.exit.A, 0.5, options);
    nozzle.PtRatio = nozzle.throat.A/nozzle.exit.A/AreaMachFunc(gam,exit.M);
    shock.M = fsolve( @(x) (((gam+1)*x^2/2)/(1 + (gam-1)*x^2/2))^(gam/(gam-1))*(((gam+1)/2)/(gam*x^2 - (gam-1)/2))^(1/(gam-1)) - nozzle.PtRatio,2,options);
    shock.A = nozzle.throat.A/AreaMachFunc(gam,shock.M);
    shock.x = fsolve( @(x) A(x) - shock.A, (nozzle.geometry.length + nozzle.geometry.xThroat)/2, options);
    shockInNozzle = true;
    fprintf('WARNING: Shock in nozzle not enabled for non-ideal nozzle.\n');
elseif (pressureRatio < exit.critical.PtRatioSupersonic - deltaPtRatio)
    %fprintf('ideal: Overexpanded flow\n');
    nozzle.status = 'ideal: overexpanded';
elseif (pressureRatio < exit.critical.PtRatioSupersonic + deltaPtRatio)
    %fprintf('Approximately fully expanded flow\n');
    nozzle.status = 'ideal: fully expanded';
else
    %fprintf('ideal: Underexpanded flow\n');
    nozzle.status = 'ideal: underexpanded';
end

% ======================= NOZZLE FRICTION & HEAT =========================

% Make a first guess at stagnation temperature and Cf along nozzle
dTstagdx = @(x) 0*x;
Tstag = @(x) inlet.Tstag;
Cf = @(x) 0.002;

xPositionOld = [0; nozzle.geometry.length];
flow.Tstag = [inlet.Tstag; inlet.Tstag];

converged = false;
maxIterations = 10; % max number of iterations to solve for Cf and Tstag
counter = 0; % used to count number of iterations
tolerance = error.betweenIterations.exitTemp; % tolerance for percent error in exit static temperature between iterations
Texit_old = 0; % saves previous exit static temperature

while ~converged
    
    counter = counter + 1;
    
    % ========================== SOLVE 1-D E.O.M =============================
    if (strcmp(nozzle.status,'ideal: subsonic')) % Subsonic flow present in nozzle
        dM2dxSubsonic = @(x, M2) -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(nozzle.geometry.length-x)./A(nozzle.geometry.length-x) + 2*gam*M2*Cf(nozzle.geometry.length-x)./D(nozzle.geometry.length-x) + (1+gam*M2)*dTstagdx(nozzle.geometry.length-x)./(2*Tstag(nozzle.geometry.length-x)));
        
        % ODE solver options
        options.RelTol = 1e-8;
        options.AbsTol = 1e-8;
        
        % Solve using 4th-order Runge-KuTstaga method
        [xPosition,M2] = ode45(dM2dxSubsonic,[0 nozzle.geometry.length],exit.M^2,options);

        xPosition = -flipud(xPosition) + nozzle.geometry.length;
        M2 = flipud(M2);

        PstagExit = inlet.Pstag*(nozzle.inlet.A/nozzle.exit.A).*(AreaMachFunc(gam,sqrt(M2(1)))./AreaMachFunc(gam,exit.M)).*sqrt(Tstag(nozzle.geometry.length)/inlet.Tstag); % stagnation pressure from mass conservation
        Pexit = PstagExit/(1 + (gam-1)*M2(end)/2).^(gam/(gam-1));
        
        fprintf('! Subsonic calculations for non-ideal nozzle not set up correctly yet\n');
        
%         % Iterate if Pexit ~= freestream.P
%         tolerance2 = 1e-6; % percent error tolerated for convergence of subsonic flow calcs
%         counter2 = 0; % used in number of subsonic iterations
%         while ( abs(Pexit - freestream.P)/freestream.P > tolerance2)            
%             
%             % Perform a Newton iteration
%             Pfunc = (Pexit - freestream.P)/freestream.P;
%             delta = 1e-12;
%             
%             % Calculate derivative
%             [~,M2B] = ode45(dM2dxSubsonic,[0 nozzle.xExit],(exit.M+delta)^2,options);
%             M2B = flipud(M2B);
%             PstagExitB = inlet.Pstag*(nozzle.inlet.A/nozzle.exit.A).*(AreaMachFunc(gam,sqrt(M2B(1)))./AreaMachFunc(gam,(exit.M+delta))).*sqrt(Tstag(nozzle.xExit)/inlet.Tstag); % stagnation pressure from mass conservation
%             PexitB = PstagExitB/(1 + (gam-1)*M2B(end)/2).^(gam/(gam-1));
%             Pderiv = (PexitB - Pexit)/delta/freestream.P;
%             
%             % Calculate new exit Mach number
%             exit.M = exit.M - Pfunc/Pderiv;
%             
%             % Solve using 4th-order Runge-KuTstaga method
%             [xPosition,M2] = ode45(dM2dxSubsonic,[0 nozzle.xExit],exit.M^2,options);
% 
%             xPosition = -flipud(xPosition) + nozzle.xExit;
%             M2 = flipud(M2);
% 
%             PstagExit = inlet.Pstag*(nozzle.inlet.A/nozzle.exit.A).*(AreaMachFunc(gam,sqrt(M2(1)))./AreaMachFunc(gam,exit.M)).*sqrt(Tstag(nozzle.xExit)/inlet.Tstag); % stagnation pressure from mass conservation
%             Pexit = PstagExit/(1 + (gam-1)*M2(end)/2).^(gam/(gam-1));
%             
%             counter2 = counter2 + 1;
%             
%             % Terminate if more than 3 iterations used
%             if(counter2 >= 3)
%                 break;
%             end
%             
%         end
        
    else % Supersonic flow present in nozzle
        
        % Estimate dMdx in order to better smooth the transition b/w
        % subsonic and supersonic flow.
        dMdxCoeff = -dAdx(nozzle.geometry.xThroat)./A(nozzle.geometry.xThroat) + 2*gam*Cf(nozzle.geometry.xThroat)./D(nozzle.geometry.xThroat) + (1+gam)*dTstagdx(nozzle.geometry.xThroat)./(2*Tstag(nozzle.geometry.xThroat));
        if(dMdxCoeff < 0)
            dMdxCoeff = -dMdxCoeff;
        end
        
        % Solve for location where dMdx = 0 (location of apparent throat)
        dMdxCoeffFunc = @(x) -dAdx(x)./A(x) + 2*gam*Cf(x)./D(x) + (1+gam)*dTstagdx(x)./(2*Tstag(x));
        options2 = optimset('TolFun',error.solver.apparentThroatLocation);
        dMdxCoeffFunc(nozzle.geometry.xThroat);
        nozzle.geometry.xApparentThroat = fzero(dMdxCoeffFunc,nozzle.geometry.xThroat,options2);
        
        % Split problem into before and after nozzle throat, solve for d(M^2)/dx
        dM2dxPost = @(x, M2) (2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(x+nozzle.geometry.xApparentThroat)./A(x+nozzle.geometry.xApparentThroat) + 2*gam*M2*Cf(x+nozzle.geometry.xApparentThroat)./D(x+nozzle.geometry.xApparentThroat) + (1+gam*M2)*dTstagdx(x+nozzle.geometry.xApparentThroat)./(2*Tstag(x+nozzle.geometry.xApparentThroat)));
        dM2dxPrior = @(x, M2) -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(nozzle.geometry.xApparentThroat-x)./A(nozzle.geometry.xApparentThroat-x) + 2*gam*M2*Cf(nozzle.geometry.xApparentThroat-x)./D(nozzle.geometry.xApparentThroat-x) + (1+gam*M2)*dTstagdx(nozzle.geometry.xApparentThroat-x)./(2*Tstag(nozzle.geometry.xApparentThroat-x)));
        
        % Make a heuristic estimate of dMdx at apparent throat
        dMdx = 600*dMdxCoeff/error.dMdxDenominator; % 600 corresponds dMdx for linear interpolation between M = 0.999 and M = 1.001
        UpperM = 1.001; % start integration at this Mach number for aft portion of nozzle
        LowerM = 0.999; % start integration at this Mach number for fore portion of nozzle
        %fprintf('dMdx: %f\n',dMdx);
        
        % ODE solver options
        options.RelTol = error.solver.M2relative;
        options.AbsTol = error.solver.M2absolute;
        % Solve using 4th-order Runge-Kutta method
        [xPositionPost,M2Post] = ode45(dM2dxPost,[(UpperM-1)/dMdx nozzle.geometry.length - nozzle.geometry.xApparentThroat],UpperM.^2,options);
        [xPositionPrior,M2Prior] = ode45(dM2dxPrior,[(1-LowerM)/dMdx nozzle.geometry.xApparentThroat],LowerM.^2,options);
        
        % Combine both parts of problem
        M2 = [flipud(M2Prior); M2Post]; % contains Mach^2
        xPosition = [-flipud(xPositionPrior)+nozzle.geometry.xApparentThroat; xPositionPost + nozzle.geometry.xApparentThroat];
    
    end
    
    % Throw exception if M^2 is negative for whatever reason
    if (sum(M2<0) > 0) % i.e. if any M2 < 0, implying Mach number is imaginary
        error('! Imaginary Mach number calculated')
    end
    
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
    %Tstag = @(x) interp1(xPosition,flow.Tstag,x,'linear');
    dTstagdxVal = (freestream.T - flow.Tstag)*4./(Cp(flow.T).*flow.density.*flow.U.*D(xPosition).*(1./flow.hf + t(xPosition)/nozzle.wall.k + 1/nozzle.hInf));
    dTstagdx = @(x) interpLinear(xPosition,dTstagdxVal,x);
    %dTstagdx = @(x) interp1(xPosition,dTstagdxVal,x,'linear');
    
    % Estimate interior wall temperature
    Qw = Cp(flow.T).*flow.density.*flow.U.*D(xPosition).*dTstagdxVal/4;
    nozzle.Tw = flow.Tstag + Qw./flow.hf; % wall temperature
    nozzle.wallRecoveryFactor = (nozzle.Tw./flow.T - 1)./((gam-1)*flow.M.^2/2);
    
    % Estimate exterior wall temperature
    nozzle.Text = ones(length(xPosition),1)*freestream.T - Qw./nozzle.hInf; % nozzle exterior wall temp.
    
    % Redefine friction coefficient distribution (Sommer & Short's method)
    TPrimeRatio = 1 + 0.035*flow.M.^2 + 0.45*(nozzle.Tw./flow.T -1);
    RePrimeRatio = 1./(TPrimeRatio.*(TPrimeRatio).^1.5.*(1 + 110.4./flow.T)./(TPrimeRatio + 110.4./flow.T));
    CfIncomp = 0.074./flow.Re.^0.2;
    flow.Cf = CfIncomp./TPrimeRatio./RePrimeRatio.^0.2;
    Cf = @(x) interpLinear(xPosition,flow.Cf,x);
    %Cf = @(x) interp1(xPosition,flow.Cf,x,'linear');
    
    % Save old solution x position for next iteration
    xPositionOld = xPosition;
    
    % =========== ESTIMATE ACCURACY OF CRITICAL FLOW LOCATION ================
    % This is still yet to be implemented, if it is at all important.
    
    fprintf('%i ',counter);
    
    if counter >= maxIterations
        fprintf('! Max iterations reached for non-ideal nozzle heat transfer & friction.\n');
        break;
    end
    
    % Check tolerance on static temperature at nozzle exit
    percentError = abs(T(nozzle.geometry.length)-Texit_old)/T(nozzle.geometry.length);
    Texit_old = T(nozzle.geometry.length);
    if(percentError < tolerance)
        converged = true;
        fprintf('iter to converge non-ideal nozzle heat xfer & friction calcs\n');
    end
    
end

% Assign flow data to nozzle 
nozzle.flow = flow;

nozzle.PstagRatio = flow.Pstag(end)/flow.Pstag(1);
nozzle.exit.Pstag = flow.Pstag(end);
nozzle.TstagRatio = flow.Tstag(end)/flow.Tstag(1);
nozzle.exit.Tstag = flow.Tstag(end);

% ========================== CALC STRESSES ===============================

nozzle.hoopStress = flow.P.*D(xPosition)./(2*t(xPosition));
nozzle.thermalHoopStress = 0.5*(nozzle.Tw-nozzle.Text)*nozzle.wall.coeffThermalExpansion*nozzle.wall.E/(1-nozzle.wall.poissonRatio);
%nozzle.thermalLongitudinalStress = 0.5*tempDiff*nozzle.wall.coeffThermalExpansion*nozzle.wall.E/(1-nozzle.wall.poissonRatio);

nozzle.maxStress = nozzle.hoopStress + nozzle.thermalHoopStress;

% ========================== CALC GEOMETRY ===============================

nozzle.xPosition = xPosition;
nozzle.geometry.A = A(xPosition);
nozzle.geometry.dAdx = dAdx(xPosition);
nozzle.geometry.D = D(xPosition);
nozzle.wall.t = t(xPosition);

% =================== CALC NOZZLE MATERIAL VOLUME ========================

% Volume calculation only works for spline parameterized nozzle geometry
% and piecewise-linear parameterized nozzle wall thickness
nozzle.geometry.volume = wallVolume(pp,nozzle.wall);

end


