% clc;
clear all; close all;

% === Dakota interface script for nozzleIdeal.m and nozzleNonIdeal.m ======
%
% Rick Fenrich 7/29/15 modified 10/16/15
% Modified by Jason Monschke 10/27/15

% ========================== INPUT PARAMETERS =============================

altitude = NaN;
mach = NaN;

% ------------------- SET ERROR TOLERANCE RANGES -------------------------
% Set error tolerances for various iterations and solvers
error.betweenIterations.exitTemp = 1e-8;
error.solver.apparentThroatLocation = 1e-6;
error.solver.M2relative = 1e-10;
error.solver.M2absolute = 1e-10;
error.dMdxDenominator = 4; % this is not an error tolerance, rather it is used to set the slope of dMdx in the transonic regime

% Read parameters from Dakota generated parameter file:
path = textscan(pwd,'%s','delimiter','.');
mytag = NaN;
if(numel(path{:})>1)
    mytag_temp = path{:}(end);
    mytag = str2num(mytag_temp{:});
end

if(isnan(mytag))
    parameters_file = 'params.in';
    results_file = 'results.out';
else
    parameters_file = sprintf('params.in.%d',mytag);
    results_file = sprintf('results.out.%d',mytag);
end

fid = fopen(parameters_file);
tline = fgets(fid);
a = sscanf(tline, ' %d %s\n');
num_vars = a(1);
vars_text = char(a(2:end).');

for ii = 1:num_vars
    tline = fgets(fid);
    a = sscanf(tline, ' %f %s\n');
    var_i = a(1);
    label_i = char(a(2:end).');
    
    if(strcmp(label_i,'alt'))
        altitude = var_i;
    elseif(strcmp(label_i,'mach'))
        mach = var_i;
    elseif(strcmp(label_i,'hInf'))
        nozzle.hInf = var_i;
    elseif(strcmp(label_i,'kWall'))
        nozzle.wall.k = var_i;
    elseif(strcmp(label_i,'TstagIn'))
        nozzle.inlet.Tstag = var_i;
    elseif(strcmp(label_i,'PstagIn'))
        nozzle.inlet.Pstag = var_i;
    elseif(strcmp(label_i,'t1'))
        t1 = var_i;
    elseif(strcmp(label_i,'t2'))
        t2 = var_i;
    elseif(strcmp(label_i,'t3'))
        t3 = var_i;
    elseif(strcmp(label_i,'rThroat'))
        rThroat = var_i;
    elseif(strcmp(label_i,'rExit'))
        rExit = var_i;
    else
        disp('Unknown variable!');
    end
end
fclose(fid);


if(isnan(altitude))
    altitude = 0;
end

if(isnan(mach))
    mach = 0;
end


atm = StndAtm(altitude*0.3048,'SI');
freestream.P = atm.P; % Pa, atmospheric pressure
freestream.T = atm.T; % K, atmospheric temperature

gam = 1.4; % ratio of specific heats
Cp = 1006; % J/kg*K, specific heat for dry air
R = 287.06; % J/kg*K, specific gas constant for dry air

% -------------------------- USEFUL FUNCTIONS ----------------------------

AreaMachFunc = @(g,M) ((g+1)/2)^((g+1)/(2*(g-1)))*M./(1+(g-1)*M.^2/2).^((g+1)/(2*(g-1)));
massFlowRate = @(Pstag,Area,Tstag,Mach) (gam/((gam+1)/2)^((gam+1)/(2*(gam-1))))*Pstag*Area*AreaMachFunc(gam,Mach)/sqrt(gam*R*Tstag);

freestream.M = mach;
freestream.U = freestream.M*sqrt(gam*R*freestream.T);


% ---------------------- SET NOZZLE GEOMETRY -----------------------------
nozzle.geometry.shape = 'spline'; % options include 'linear' and 'spline'
nozzle.geometry.xThroat = 0.33;
nozzle.geometry.length = 1;
rInlet = 0.3255;
% rThroat = 0.278296628003098;
% rExit = 0.329285010926785;

nozzle.inlet.D = rInlet*2;
nozzle.geometry.Ainlet2Athroat = rInlet^2/rThroat^2;
nozzle.geometry.Aexit2Athroat = rExit^2/rThroat^2;

nozzle.geometry.spline.seed = [0    rInlet;
                        nozzle.geometry.xThroat rThroat;
                        nozzle.geometry.length    rExit];
nozzle.geometry.spline.slopes = [0, 0];

% -------------------- SET NOZZLE WALL GEOMETRY --------------------------
nozzle.wall.shape = 'piecewise-linear';
% nozzle.wall.seed = an array of the form [x; y] where [x,y]
% denote the location of the control points with the origin being at
% the center of the inlet area
% nozzle.wall.breaks = vector giving location of breaks in piecewise
% function

nozzle.wall.seed = [0, t1; 
                    nozzle.geometry.xThroat, t2; 
                    nozzle.geometry.length, t3];
nozzle.wall.breaks = [0;
                      nozzle.geometry.xThroat;
                      nozzle.geometry.length];

% --------------------- FLUID INPUT PROPERTIES ---------------------------
fluid.gam = 1.4; % ratio of specific heats
fluid.R = 287.06; % J/kg-K, specific gas constant

% ====================== RUN NOZZLE CALCULATIONS =========================

[ nozzle ] = nozzleNonIdeal( fluid, freestream, nozzle, error );

% ============================ OUTPUT DATA ===============================

% Nozzle exit parameters
nozzle.exit.Pstag = nozzle.inlet.Pstag*nozzle.PstagRatio;
nozzle.exit.M = nozzle.flow.M(end);
nozzle.exit.T = nozzle.flow.T(end);
nozzle.exit.U = nozzle.exit.M*sqrt(gam*R*nozzle.exit.T);
nozzle.exit.P = nozzle.exit.Pstag/(1 + (gam-1)*nozzle.exit.M^2/2)^(gam/(gam-1));

nozzle.massFlowRate = massFlowRate(nozzle.inlet.Pstag,nozzle.inlet.A,nozzle.inlet.Tstag,nozzle.flow.M(1));

thrust = nozzle.massFlowRate*(nozzle.exit.U - freestream.U) + (nozzle.exit.P - freestream.P)*nozzle.exit.A;

fid = fopen(results_file,'w');
fprintf(fid,'%.16e thrust\n',thrust);
fprintf(fid,'%.16e maxStress\n',max(nozzle.hoopStress));
fprintf(fid,'%.16e maxTemp\n',max(nozzle.Tw));
fprintf(fid,'%.16e volume\n',nozzle.geometry.volume);
fclose(fid);


