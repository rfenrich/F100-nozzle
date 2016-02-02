function [ nozzle ] = nozzleIdeal( fluid, freestream, nozzle, error)
% Solve for flow along length of ideal nozzle given geometry, inlet
% stagnation temperature and pressure, and freestream temperature and
% pressure. Returns M, density, pressure P, temperature T, stagnation 
% temp. Tstag, stagnation pressure Pstag, velocity U, Re, and hoopStress
% along length of nozzle.
%
% INPUTS:
% fluid = structure with fields: gam (ratio of specific heats) and R
%         (specific ideal gas constant)
% freestream = structure with fields: T and P
% nozzle = structure with fields: geometry.Ainlet2Athroat, 
%          geometry.Aexit2Athroat, geometry.length, geometry.shape,
%          geometry.xThroat,  geometry.xExit, inlet.Tstag, inlet.Pstag,
%          inlet.D
% hInf = generalized heat transfer coeff. from outside nozzle wall to
%        freestream (units W/m^2-K)
% error = structure defining error tolerances for iterations and solvers
%         included the following fields: error.solver.M2relative, 
%         error.solver.M2absolute
%
% OUTPUTS:
% nozzle = modified input structure with additional fields including flow
% and specific geometry
%
% Rick Fenrich 7/10/15 modified 10/16/15

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

% Calculate pressure ratio which determines state of nozzle:
pressureRatio = inlet.Pstag/freestream.P;

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
    
elseif(strcmp(nozzle.geometry.shape,'B-spline'))
    
    % Adjust nozzle throat size/location information if it has changed
    [xThroat, yThroat] = BsplineGeometry(0, 'throat', nozzle.geometry.bSpline.knots, nozzle.geometry.bSpline.coefs);
    if(xThroat ~= nozzle.geometry.xThroat)
        fprintf('throat size/location changed with spline parameterization\n');
    end
    nozzle.geometry.xThroat = xThroat;
    nozzle.throat.A = pi*yThroat^2;
    nozzle.geometry.Ainlet2Athroat = nozzle.inlet.A/nozzle.throat.A;
    nozzle.geometry.Aexit2Athroat = nozzle.exit.A/nozzle.throat.A;

    % Make necessary functions for splined nozzle shape
    A = @(x) BsplineGeometry(x,'A',nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs);
    dAdx = @(x) BsplineGeometry(x,'dAdx',nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs);
    D = @(x) BsplineGeometry(x,'D',nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs);

elseif(strcmp(nozzle.geometry.shape,'B-spline-mex'))
    
    % Adjust nozzle throat size/location information if it has changed
    [xThroat, yThroat] = BsplineGeometry(0, 'throat', nozzle.geometry.bSpline.knots, nozzle.geometry.bSpline.coefs);
    if(xThroat ~= nozzle.geometry.xThroat)
        fprintf('throat size/location changed with spline parameterization\n');
    end
    nozzle.geometry.xThroat = xThroat;
    nozzle.throat.A = pi*yThroat^2;
    nozzle.geometry.Ainlet2Athroat = nozzle.inlet.A/nozzle.throat.A;
    nozzle.geometry.Aexit2Athroat = nozzle.exit.A/nozzle.throat.A;

    % Make necessary functions for splined nozzle shape
    A = @(x) BsplineGeometryMex(x,1,nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs');
    dAdx = @(x) BsplineGeometryMex(x,2,nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs');
    D = @(x) BsplineGeometryMex(x,3,nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs');

else % if nozzle shape is not a spline
    A = @(x) nozzleGeometry(x,'A',nozzle.inlet.D,nozzle.geometry.length,nozzle.geometry.xThroat,nozzle.geometry.Ainlet2Athroat,nozzle.geometry.Aexit2Athroat,nozzle.geometry.shape);
    dAdx = @(x) nozzleGeometry(x,'dAdx',nozzle.inlet.D,nozzle.geometry.length,nozzle.geometry.xThroat,nozzle.geometry.Ainlet2Athroat,nozzle.geometry.Aexit2Athroat,nozzle.geometry.shape);
    D = @(x) nozzleGeometry(x,'D',nozzle.inlet.D,nozzle.geometry.length,nozzle.geometry.xThroat,nozzle.geometry.Ainlet2Athroat,nozzle.geometry.Aexit2Athroat,nozzle.geometry.shape);
end

% ======================= NOZZLE WALL GEOMETRY ===========================
% Only piecewise-linear geometry is enabled thus far
t = @(x) piecewiseLinearGeometry(x,'t',nozzle.wall); % m, thickness of wall

% ======================= DETERMINE NOZZLE FLOW ==========================
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
    nozzle.status = 'subsonic';
    exit.M = sqrt(2/(gam-1))*sqrt(pressureRatio^((gam-1)/gam) - 1);
    nozzle.PstagRatio = 1;
elseif (pressureRatio < exit.normalShock.PtRatio)
    %fprintf('ideal: Shock in nozzle\n');
    nozzle.status = 'shock in nozzle';
    exit.M = fsolve( @(x) ((gam+1)/2)^((gam+1)/(2*(gam-1)))*x*sqrt(1 + (gam-1)*x^2/2) - pressureRatio*nozzle.throat.A/nozzle.exit.A, 0.5, options);
    nozzle.PstagRatio = nozzle.throat.A/nozzle.exit.A/AreaMachFunc(gam,exit.M);
    shock.M = fsolve( @(x) (((gam+1)*x^2/2)/(1 + (gam-1)*x^2/2))^(gam/(gam-1))*(((gam+1)/2)/(gam*x^2 - (gam-1)/2))^(1/(gam-1)) - nozzle.PstagRatio,2,options);
    shock.A = nozzle.throat.A/AreaMachFunc(gam,shock.M);
    shock.x = fsolve( @(x) A(x) - shock.A, (nozzle.geometry.length + nozzle.geometry.xThroat)/2, options);
    shockInNozzle = true;
elseif (pressureRatio < exit.critical.PtRatioSupersonic - deltaPtRatio)
    %fprintf('ideal: Overexpanded flow\n');
    nozzle.status = 'overexpanded';
    nozzle.PstagRatio = 1;
elseif (pressureRatio < exit.critical.PtRatioSupersonic + deltaPtRatio)
    %fprintf('Approximately fully expanded flow\n');
    nozzle.status = 'fully expanded';
    nozzle.PstagRatio = 1;
else
    %fprintf('ideal: Underexpanded flow\n');
    nozzle.status = 'underexpanded';
    nozzle.PstagRatio = 1;
end

% ========================== SOLVE 1-D E.O.M =============================
% Split problem into before and after nozzle throat, solve for d(M^2)/dx
dM2dxPost = @(x, M2) (2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(x+nozzle.geometry.xThroat)./A(x+nozzle.geometry.xThroat));
dM2dxPrior = @(x, M2) -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(nozzle.geometry.xThroat-x)./A(nozzle.geometry.xThroat-x));

% ODE solver options
options.RelTol = error.solver.M2relative;
options.AbsTol = error.solver.M2absolute;
% Solve using 4th-order Runge-Kutta method
if (strcmp(nozzle.status,'subsonic')) % subsonic flow throughout nozzle
    dM2dxSubsonic = @(x, M2) -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(nozzle.geometry.length-x)./A(nozzle.geometry.length-x));
                             
    [xPosition,M2] = ode45(dM2dxSubsonic,[0 nozzle.geometry.length],exit.M^2,options);
    xPosition = -flipud(xPosition) + nozzle.geometry.length;
    M2 = flipud(M2);
elseif (shockInNozzle ~= true) % supersonic flow, no shock in nozzle
    
    UpperM = 1.000001; % start integration at this Mach number for aft portion of nozzle
    LowerM = 0.999999; % start integration at this Mach number for fore portion of nozzle
    
    % Solve using 4th-order Runge-Kutta method
    [xPositionPost,M2Post] = ode45(dM2dxPost,[0 nozzle.geometry.length-nozzle.geometry.xThroat],UpperM,options);
    [xPositionPrior,M2Prior] = ode45(dM2dxPrior,[0 nozzle.geometry.xThroat],LowerM,options);
    
    % Combine both parts of problem
    M2 = [flipud(M2Prior); M2Post]; % contains Mach^2
    xPosition = [-flipud(xPositionPrior)+nozzle.geometry.xThroat; xPositionPost + nozzle.geometry.xThroat];
else % sub and supersonic flow, shock in nozzle
    dM2dxPostShock = @(x, M2) (2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(x+shock.x)./A(x+shock.x));
    
    UpperM = 1.000001; % start integration at this Mach number for aft portion of nozzle
    LowerM = 0.999999; % start integration at this Mach number for fore portion of nozzle
    
    [xPositionPost,M2Post] = ode45(dM2dxPost,[0 shock.x - nozzle.geometry.xThroat],UpperM,options);
    MbehindShock = sqrt((1 + (gam-1)*shock.M^2/2)/(gam*shock.M^2 - (gam-1)/2));
    shock.PtRatio = nozzle.PstagRatio;
    [xPositionPostShock,M2PostShock] = ode45(dM2dxPostShock,[1e-8 nozzle.geometry.length-shock.x],MbehindShock^2,options);
    [xPositionPrior,M2Prior] = ode45(dM2dxPrior,[0 nozzle.geometry.xThroat],LowerM,options);
    
    % Combine both parts of problem
    M2 = [flipud(M2Prior); M2Post; M2PostShock]; % contains Mach^2
    xPosition = [-flipud(xPositionPrior)+nozzle.geometry.xThroat; xPositionPost + nozzle.geometry.xThroat; xPositionPostShock + shock.x];
end

% Throw exception if M^2 is negative for whatever reason
if (sum(M2<0) > 0) % i.e. if any M2 < 0, implying Mach number is imaginary
    error('! Imaginary Mach number calculated')
end

% ======================= CALC OTHER PROPERTIES ==========================

flow.M = sqrt(M2); % Mach number
flow.Tstag = inlet.Tstag*ones(length(xPosition),1);
flow.T = flow.Tstag./(1 + (gam-1)*M2/2); % static temperature from stag. temp. definition

if (shockInNozzle ~= true) % no shock in nozzle
    flow.Pstag = inlet.Pstag*ones(length(xPosition),1);
    nozzle.exit.Pstag = inlet.Pstag;
else % shock in nozzle
    flow.Pstag = [inlet.Pstag*ones(length(xPositionPrior)+length(xPositionPost),1); inlet.Pstag*shock.PtRatio*ones(length(xPositionPostShock),1)];
    nozzle.exit.Pstag = flow.Pstag(end);
end

flow.P = flow.Pstag./(1 + (gam-1)*M2/2).^(gam/(gam-1)); % static pressure from stag. press. definition
flow.density = flow.P./(R*flow.T); % density
flow.U = flow.M.*sqrt(gam*R*flow.T); % velocity
flow.Re = flow.density.*flow.U.*D(xPosition)./dynamicViscosity(flow.T);

% =================== ASSIGN FLOW DATA TO NOZZLE =========================

nozzle.flow = flow;
nozzle.PstagRatio = flow.Pstag(end)/flow.Pstag(1);
nozzle.TstagRatio = 1;

% Record nozzle exit values
nozzle.exit.M = nozzle.flow.M(end);
nozzle.exit.Pstag = flow.Pstag(end);
nozzle.exit.Tstag = inlet.Tstag;
nozzle.exit.P = nozzle.flow.P(end);
nozzle.exit.T = nozzle.flow.T(end);
nozzle.exit.U = nozzle.flow.U(end);

% ==================== CALCULATE APPROX. THRUST ==========================
% Thrust is only approximate since engine is not taken into account.

massFlowRate = @(Pstag,Area,Tstag,M) (gam/((gam+1)/2)^((gam+1)/(2*(gam-1))))*Pstag*Area*AreaMachFunc(gam,M)/sqrt(gam*R*Tstag);
nozzle.massFlowRate = massFlowRate(nozzle.inlet.Pstag,nozzle.inlet.A,nozzle.inlet.Tstag,nozzle.flow.M(1));
nozzle.approxThrust = nozzle.massFlowRate*(nozzle.exit.U - freestream.U) + (nozzle.exit.P - freestream.P)*nozzle.exit.A;

% ========================== CALC STRESSES ===============================

nozzle.hoopStress = flow.P.*D(xPosition)./(2*t(xPosition));

% ========================== CALC GEOMETRY ===============================

nozzle.xPosition = xPosition;
nozzle.geometry.A = A(xPosition);
nozzle.geometry.dAdx = dAdx(xPosition);
nozzle.geometry.D = D(xPosition);
nozzle.wall.t = t(xPosition);

% =================== CALC NOZZLE MATERIAL VOLUME ========================

% Volume calculation only works for spline parameterized nozzle geometry
% and piecewise-linear parameterized nozzle wall thickness
if(exist('pp','var')) % Exact volume for cubic spline parameterization
    nozzle.geometry.volume = wallVolume(pp,nozzle.wall);
else % Approximate volume using trapezoidal integration
    xVolume = linspace(0,nozzle.geometry.length,500)';
    volumeIntegrand = pi*D(xVolume).*t(xVolume) + pi*t(xVolume).^2;
    nozzle.geometry.volume = (xVolume(2)-xVolume(1))*trapz(volumeIntegrand);
end

end
