function [ flow, nozzle, xPosition ] = nozzleIdeal( fluid, inlet, freestream, nozzle)
% Solve for flow along length of ideal nozzle given geometry, inlet
% stagnation temperature and pressure, and freestream temperature and
% pressure. Returns M, density, pressure P, temperature T, stagnation 
% temp. Tstag, stagnation pressure Pstag, velocity U, Re, and hoopStress
% along length of nozzle.
%
% INPUTS:
% fluid = structure with fields: gam (ratio of specific heats) and R
%         (specific ideal gas constant)
% inlet = structure with fields: Tstag, Pstag, D (diameter)
% freestream = structure with fields: T and P
% nozzle = structure with fields: Ainlet2Athroat, Aexit2Athroat, length, shape,
%          xThroat,  xExit
%
% OUTPUTS:
% flow = structure with vectors of flow properties along length of nozzle
% nozzle = modified input structure with additional fields
% xPosition = vector denoting x-location of flow properties in flow struct
%
% Rick Fenrich 7/10/15

% ========================== GAS PROPERTIES ==============================
gam = fluid.gam;
R = fluid.R;

% Area-Mach function from mass 1-D mass conservation equations:
AreaMachFunc = @(g,M) ((g+1)/2)^((g+1)/(2*(g-1)))*M./(1+(g-1)*M.^2/2).^((g+1)/(2*(g-1)));
% Sutherland's law for dynamic viscosity:
dynamicViscosity = @(T) 1.716e-5*(T/273.15).^1.5*(273.15 + 110.4)./(T + 110.4); % kg/m*s

% ========================= NOZZLE PROPERTIES ============================
% Calculate nozzle inlet, throat, and exit areas if they are not given:
if(~exist('nozzle.inlet.A','var'))
    nozzle.inlet.A = pi*nozzle.inlet.D^2/4;
end
if(~exist('nozzle.throat.A','var'))
    nozzle.throat.A = nozzle.inlet.A/nozzle.Ainlet2Athroat;
end
if(~exist('nozzle.exit.A','var'))
    nozzle.exit.A = nozzle.Aexit2Athroat*nozzle.inlet.A/nozzle.Ainlet2Athroat;
end

% Calculate pressure ratio which determines state of nozzle:
pressureRatio = inlet.Pstag/freestream.P;

% ========================== NOZZLE GEOMETRY =============================

A = @(x) nozzleGeometry(x,'A',nozzle.inlet.D,nozzle.length,nozzle.xThroat,nozzle.Ainlet2Athroat,nozzle.Aexit2Athroat,nozzle.shape);
dAdx = @(x) nozzleGeometry(x,'dAdx',nozzle.inlet.D,nozzle.length,nozzle.xThroat,nozzle.Ainlet2Athroat,nozzle.Aexit2Athroat,nozzle.shape);
D = @(x) nozzleGeometry(x,'D',nozzle.inlet.D,nozzle.length,nozzle.xThroat,nozzle.Ainlet2Athroat,nozzle.Aexit2Athroat,nozzle.shape);
t = @(x) nozzleGeometry(x,'t',nozzle.inlet.D,nozzle.length,nozzle.xThroat,nozzle.Ainlet2Athroat,nozzle.Aexit2Athroat,nozzle.shape); % m, thickness of wall

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
    shock.x = fsolve( @(x) A(x) - shock.A, (nozzle.xExit + nozzle.xThroat)/2, options);
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
dM2dxPost = @(x, M2) (2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(x+nozzle.xThroat)./A(x+nozzle.xThroat));
dM2dxPrior = @(x, M2) -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(nozzle.xThroat-x)./A(nozzle.xThroat-x));

% ODE solver options
options.RelTol = 1e-8;
options.AbsTol = 1e-8;
% Solve using 4th-order Runge-Kutta method
if (strcmp(nozzle.status,'subsonic')) % subsonic flow throughout nozzle
    dM2dxSubsonic = @(x, M2) -(2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(nozzle.xExit-x)./A(nozzle.xExit-x));
                             
    [xPosition,M2] = ode45(dM2dxSubsonic,[0 nozzle.xExit],exit.M^2,options);
    xPosition = -flipud(xPosition) + nozzle.xExit;
    M2 = flipud(M2);
elseif (shockInNozzle ~= true) % supersonic flow, no shock in nozzle
    [xPositionPost,M2Post] = ode45(dM2dxPost,[1e-6 nozzle.xExit-nozzle.xThroat],1.005,options);
    [xPositionPrior,M2Prior] = ode45(dM2dxPrior,[1e-6 nozzle.xThroat],0.999,options);
    
    % Combine both parts of problem
    M2 = [flipud(M2Prior); M2Post]; % contains Mach^2
    xPosition = [-flipud(xPositionPrior)+nozzle.xThroat; xPositionPost + nozzle.xThroat];
else % sub and supersonic flow, shock in nozzle
    dM2dxPostShock = @(x, M2) (2*M2*(1+(gam-1)*M2/2)/(1-M2))*(-dAdx(x+shock.x)./A(x+shock.x));
    
    [xPositionPost,M2Post] = ode45(dM2dxPost,[1e-6 shock.x-nozzle.xThroat],1.001,options);
    MbehindShock = sqrt((1 + (gam-1)*shock.M^2/2)/(gam*shock.M^2 - (gam-1)/2));
    shock.PtRatio = nozzle.PstagRatio;
    [xPositionPostShock,M2PostShock] = ode45(dM2dxPostShock,[1e-6 nozzle.xExit-shock.x],MbehindShock^2,options);
    [xPositionPrior,M2Prior] = ode45(dM2dxPrior,[1e-6 nozzle.xThroat],0.999,options);
    
    % Combine both parts of problem
    M2 = [flipud(M2Prior); M2Post; M2PostShock]; % contains Mach^2
    xPosition = [-flipud(xPositionPrior)+nozzle.xThroat; xPositionPost + nozzle.xThroat; xPositionPostShock + shock.x];
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

nozzle.TstagRatio = 1;
nozzle.exit.Tstag = inlet.Tstag;

% ========================== CALC STRESSES ===============================
nozzle.hoopStress = flow.P.*D(xPosition)./(2*t(xPosition));

end
