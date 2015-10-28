clc; clear all; close all;


path = './output/';
filename = 'nozzleDUU.out';

material_density = 2700; % kg/m^3

altitude = 0;
mach = 0;
nozzle.hInf = 500; % W/m^2-K, heat transfer coefficient from external nozzle wall to environment
nozzle.wall.k = 30; % W/m*K, thermal conductivity of nozzle wall
nozzle.inlet.Tstag = 888.3658;
nozzle.inlet.Pstag = 3.0550e5;

% ------------------- SET ERROR TOLERANCE RANGES -------------------------
% Set error tolerances for various iterations and solvers
error.betweenIterations.exitTemp = 1e-8;
error.solver.apparentThroatLocation = 1e-6;
error.solver.M2relative = 1e-10;
error.solver.M2absolute = 1e-10;
error.dMdxDenominator = 4; % this is not an error tolerance, rather it is used to set the slope of dMdx in the transonic regime


fileID = fopen(sprintf('%s%s',path,filename),'r');
% Search for best parameters:
while (~feof(fileID))
    line = fgetl(fileID);
    if(~isempty(findstr(line,'<<<<< Best parameters')))
        % found it!
        break;
    end
end

Nparam = 5;
for ii = 1:Nparam
    line = fgetl(fileID);
    temp = textscan(line,'%s');
    val = str2num(char(temp{1}(1)));
    if(~isempty(findstr(line,'t1')))
        t1 = val;
    elseif(~isempty(findstr(line,'t2')))
        t2 = val;
    elseif(~isempty(findstr(line,'t3')))
        t3 = val;
    elseif(~isempty(findstr(line,'rThroat')))
        rThroat = val;
    elseif(~isempty(findstr(line,'rExit')))
        rExit = val;
    end
end
fclose(fileID);

addpath('../');
atm = StndAtm(altitude*0.3048,'SI');
rmpath('../');
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

addpath('../');
[ nozzle ] = nozzleNonIdeal( fluid, freestream, nozzle, error );
rmpath('../');

% ============================ OUTPUT DATA ===============================

% Nozzle exit parameters
nozzle.exit.Pstag = nozzle.inlet.Pstag*nozzle.PstagRatio;
nozzle.exit.M = nozzle.flow.M(end);
nozzle.exit.T = nozzle.flow.T(end);
nozzle.exit.U = nozzle.exit.M*sqrt(gam*R*nozzle.exit.T);
nozzle.exit.P = nozzle.exit.Pstag/(1 + (gam-1)*nozzle.exit.M^2/2)^(gam/(gam-1));

nozzle.massFlowRate = massFlowRate(nozzle.inlet.Pstag,nozzle.inlet.A,nozzle.inlet.Tstag,nozzle.flow.M(1));

thrust = nozzle.massFlowRate*(nozzle.exit.U - freestream.U) + (nozzle.exit.P - freestream.P)*nozzle.exit.A;


%%

fprintf('Nozzle weight = %f kg or %f lbs\n',nozzle.geometry.volume*material_density, nozzle.geometry.volume*material_density*2.20462);
fprintf('Max wall thickness = %f mm\n',max(nozzle.wall.t)*1000);
fprintf('Min wall thickness = %f mm\n',min(nozzle.wall.t)*1000);


%%
% Plot nozzle shape:

fill_between_lines = @(X,Y1,Y2,C) fill( [X' fliplr(X')],  [Y1' fliplr(Y2')], C,'EdgeColor','none');

% figure;
% hold on
% plot(nozzle.xPosition,nozzle.geometry.D/2);
% plot(nozzle.xPosition,nozzle.geometry.D/2 + nozzle.wall.t);
% plot(nozzle.xPosition,-nozzle.geometry.D/2);
% plot(nozzle.xPosition,-nozzle.geometry.D/2 - nozzle.wall.t);
% plot(nozzle.xPosition,zeros(size(nozzle.xPosition)),'k--');
% xlabel('Axial position (m)');
% title('Geometry');
% axis equal;

figure;
hold on;
fill_between_lines(nozzle.xPosition,nozzle.geometry.D/2,nozzle.geometry.D/2 + nozzle.wall.t,[0 0 1]);
fill_between_lines(nozzle.xPosition,-nozzle.geometry.D/2,-nozzle.geometry.D/2 - nozzle.wall.t,[0 0 1]);
plot(nozzle.xPosition,zeros(size(nozzle.xPosition)),'k--');
xlabel('Axial position (m)');
% title('Geometry');
axis equal;
% matlab2tikz('plots/nozzleShape.tex','standalone',true);


