% ========= Driver script for nozzleIdeal.m and nozzleNonIdeal.m =========
%
% Provides inputs to the nozzleIdeal and nozzleNonIdeal functions, runs
% them, and produces plots comparing ideal and non-ideal nozzle behavior,
% as well as non-ideal nozzle heating and friction characteristics.
%
% Rick Fenrich 7/29/15 modified 10/16/15

% ========================== INPUT PARAMETERS ============================

mission = 1; % mission number for which certain input parameters are defined below
nozzle.hInf = 500; % W/m^2-K, heat transfer coefficient from external nozzle wall to environment
nozzle.wall.k = 30; % W/m*K, thermal conductivity of nozzle wall

% Define input parameters that will change based on flight regime:
if(mission == 1) % static sea-level thrust case
    altitude = 0;
    mach = 0;
    nozzle.inlet.Tstag = 888.3658;
    nozzle.inlet.Pstag = 3.0550e5;
elseif(mission == 2) % intermediate case
    altitude = 15000;
    mach = 0.5;
    nozzle.inlet.Tstag = 942.9857;
    nozzle.inlet.Pstag = 2.3227e5;
elseif(mission == 3) % high speed, high altitude case
    altitude = 35000;
    mach = 0.9;
    nozzle.inlet.Tstag = 1021.5;
    nozzle.inlet.Pstag = 1.44925e5;
elseif(mission == 4) % case with shock in nozzle
    altitude = 0;
    mach = 0.0;
    nozzle.inlet.Tstag = 900;
    nozzle.inlet.Pstag = 1.3e5;
elseif(mission == 5) % subsonic flow in nozzle
    altitude = 0;
    mach = 0.0;
    nozzle.inlet.Tstag = 900;
    nozzle.inlet.Pstag = 1.1e5;    
end

% ------------------- SET ERROR TOLERANCE RANGES -------------------------
% Set error tolerances for various iterations and solvers
error.betweenIterations.exitTemp = 1e-8;
error.solver.apparentThroatLocation = 1e-6;
error.solver.M2relative = 1e-10;
error.solver.M2absolute = 1e-10;
error.dMdxDenominator = 4; % this is not an error tolerance, rather it is used to set the slope of dMdx in the transonic regime

% ---------------------- SET NOZZLE GEOMETRY -----------------------------
nozzle.geometry.shape = 'spline'; % options include 'linear' and 'spline'
nozzle.inlet.D = 0.651; % m 
nozzle.geometry.Ainlet2Athroat = 1.368; % area ratio of inlet to throat
nozzle.geometry.Aexit2Athroat = 1.4; % area ratio of exit to throat
nozzle.geometry.length = 1; % m
nozzle.geometry.xThroat = 0.33; % m, location of throat from inlet
if(strcmp(nozzle.geometry.shape,'spline'))
    % To parameterize using a spline, the following must be provided:
    % nozzle.spline.seed = either a shape already defined in the
    % nozzleGeometry.m file or an array of the form [x; y] where [x,y]
    % denote the location of the control points with the origin being at
    % the center of the inlet area
    % nozzle.spline.nControlPoints = number of control points
    % nozzle.spline.breaks = a vector giving the x-location of
    %                                     control points
    % nozzle.spline.slopes = 1x2 array; 1st argument is slope of inlet,
    % 2nd argument is slope of outlet
    nozzle.geometry.spline.seed = ...%'linear';
                                [0, 0.3255; 
                                 0.33, 0.2783; 
                                 1, 0.3293];
    % Recalculate area ratios/inlet geometry if necessary from the spline seed:
    % nozzle.geometry.Ainlet2Athroat = (0.3255)^2/(0.2783)^2; % area ratio of inlet to throat
    % nozzle.geometry.Aexit2Athroat = (0.3293)^2/(0.2783)^2; % area ratio of exit to throat
    % nozzle.inlet.D = (0.3255)*2;
    nozzle.geometry.spline.breaks = ...
            [0;
            nozzle.geometry.xThroat;
            nozzle.geometry.length];
    nozzle.geometry.spline.slopes = [0, 0];
end

% -------------------- SET NOZZLE WALL GEOMETRY --------------------------
nozzle.wall.shape = 'piecewise-linear';
% nozzle.wall.seed = an array of the form [x; y] where [x,y]
% denote the location of the control points with the origin being at
% the center of the inlet area
% nozzle.wall.breaks = vector giving location of breaks in piecewise
% function
nozzle.wall.seed = [0, 0.01; 
                    nozzle.geometry.xThroat, 0.01; 
                    nozzle.geometry.length, 0.01];
nozzle.wall.breaks = [0;
                      nozzle.geometry.xThroat;
                      nozzle.geometry.length];

% --------------------- FLUID INPUT PROPERTIES ---------------------------
fluid.gam = 1.4; % ratio of specific heats
fluid.R = 287.06; % J/kg-K, specific gas constant

% ----------------- CALCULATE FREESTREAM PROPERTIES ----------------------
atm = StndAtm(altitude*0.3048,'SI'); % obtain standard atmosphere characteristics
freestream.P = atm.P; % Pa, atmospheric pressure
freestream.T = atm.T; % K, atmospheric temperature

% ====================== RUN NOZZLE CALCULATIONS =========================

[ nozzleI ] = nozzleIdeal( fluid, freestream, nozzle, error );

[ nozzle ] = nozzleNonIdeal( fluid, freestream, nozzle, error );

% ============================ OUTPUT DATA ===============================

fprintf('Pstag ratio: %f\n',nozzle.flow.Pstag(end)/nozzle.flow.Pstag(1));
fprintf('Tstag ratio: %f\n',nozzle.flow.Tstag(end)/nozzle.flow.Tstag(1));

% ============================== PLOTTING ================================
% Comment or uncomment the plotting below to make plots

% Manage colors and linestyles for plots:
set(groot,'defaultAxesColorOrder',[0 0 0; 1 0 0; 1 0.25 0; 0 1 0; 0 0 1])
%set(groot,'defaultAxesLineStyleOrder','-|--|:');

% -------------------- PLOT VARIOUS DATA IN 1 FIGURE ---------------------
figure
subplot(2,3,1); hold on
plot(nozzle.xPosition,nozzle.geometry.D/2)
plot(nozzle.xPosition,nozzle.geometry.D/2+nozzle.wall.t);
plot(nozzle.xPosition,-nozzle.geometry.D/2);
plot(nozzle.xPosition,-nozzle.geometry.D/2-nozzle.wall.t);
title('Geometry')
axis equal

subplot(2,3,2); hold on
plot(nozzleI.xPosition,nozzleI.flow.Re)
plot(nozzle.xPosition,nozzle.flow.Re)
title('Re')
legend('ideal');

subplot(2,3,3); hold on
plot(nozzleI.xPosition,nozzleI.flow.M)
plot(nozzle.xPosition,nozzle.flow.M)
title('Mach Number')
legend('ideal');

subplot(2,3,4); hold on
plot(nozzleI.xPosition,nozzleI.flow.P)
plot(nozzle.xPosition,nozzle.flow.P)
plot(nozzleI.xPosition,nozzleI.flow.Pstag)
plot(nozzle.xPosition,nozzle.flow.Pstag)
plot(nozzle.xPosition,freestream.P*ones(length(nozzle.xPosition),1))
title('Pressure (Pa)')
legend('Static (ideal)','Static','Stag (ideal)','Stag','\infty');

subplot(2,3,5); hold on
plot(nozzleI.xPosition,nozzleI.flow.T)
plot(nozzle.xPosition,nozzle.flow.T)
plot(nozzleI.xPosition,nozzleI.flow.Tstag)
plot(nozzle.xPosition,nozzle.flow.Tstag)
plot(nozzle.xPosition,freestream.T*ones(length(nozzle.xPosition),1),'k-')
title('Temperature (K)')
legend('Static (ideal)','Static','Stag (ideal)','Stag','\infty');

subplot(2,3,6); hold on
plot(nozzleI.xPosition,nozzleI.flow.density)
plot(nozzle.xPosition,nozzle.flow.density)
title('Density (kg/m^3)')
legend('ideal');

% % ------------------------ PLOT NOZZLE GEOMETRY --------------------------
% figure; hold on
% plot(nozzle.xPosition,nozzle.geometry.D/2)
% plot(nozzle.xPosition,-nozzle.geometry.D/2)
% xlabel('Axial position (m)')
% title('Geometry')
% axis equal
% 
% % ------------------------ PLOT REYNOLDS NUMBER --------------------------
% figure; hold on
% plot(nozzleI.xPosition,nozzleI.flow.Re)
% plot(nozzle.xPosition,nozzle.flow.Re)
% xlabel('Axial position (m)')
% title('Re')
% legend2 = legend('ideal');
% 
% % ------------------------------ PLOT MACH -------------------------------
% figure; hold on
% plot(nozzleI.xPosition,nozzleI.flow.M)
% plot(nozzle.xPosition,nozzle.flow.M)
% xlabel('Axial position (m)')
% title('Mach Number')
% legend3 = legend('ideal');
% 
% % --------------------------- PLOT PRESSURES -----------------------------
% figure; hold on
% plot(nozzleI.xPosition,nozzleI.flow.P)
% plot(nozzle.xPosition,nozzle.flow.P)
% plot(nozzleI.xPosition,nozzleI.flow.Pstag)
% plot(nozzle.xPosition,nozzle.flow.Pstag)
% plot(nozzle.xPosition,freestream.P*ones(length(nozzle.xPosition),1))
% xlabel('Axial position (m)')
% title('Pressure (Pa)')
% legend4 = legend('Static (ideal)','Static','Stag (ideal)','Stag','\infty','Location','EastOutside');
% 
% % ------------------------- PLOT TEMPERATURES ----------------------------
% figure; hold on
% hold on
% plot(nozzleI.xPosition,nozzleI.flow.T)
% plot(nozzle.xPosition,nozzle.flow.T)
% plot(nozzleI.xPosition,nozzleI.flow.Tstag)
% plot(nozzle.xPosition,nozzle.flow.Tstag)
% plot(nozzle.xPosition,freestream.T*ones(length(nozzle.xPosition),1),'k-')
% xlabel('Axial position (m)','FontName','CMU Serif','FontSize',14)
% title('Temperature (K)','FontName','CMU Serif','FontSize',16)
% set(gca,'FontName','CMU Serif','FontSize',14)
% legend5 = legend('Static (ideal)','Static','Stag (ideal)','Stag','\infty','Location','EastOutside');
% set(legend5,'FontName','CMU Serif','FontSize',12)
% 
% % ---------------------------- PLOT DENSITY ------------------------------
% figure; hold on
% plot(nozzleI.xPosition,nozzleI.flow.density)
% plot(nozzle.xPosition,nozzle.flow.density)
% xlabel('Axial position (m)')
% title('Density (kg/m^3)')
% legend6 = legend('ideal');
% 
% % --------------------------- PLOT h_f & C_f -----------------------------
% figure
% subplot(1,2,1); hold on
% plot(nozzle.xPosition,nozzle.flow.hf)
% xlabel('Axial position (m)')
% title('Wall convection coefficient h_f (W/m^2-K)')
% subplot(1,2,2); hold on
% plot(nozzle.xPosition,nozzle.flow.Cf)
% xlabel('Axial position (m)')
% title('Friction coefficient C_f')
% 
% % --------------------- PLOT TEMPERATURE PROFILES ------------------------
% figure
% hold on
% plot(nozzle.xPosition,nozzle.flow.T)
% plot(nozzle.xPosition,nozzle.flow.Tstag)
% plot(nozzle.xPosition,nozzle.Tw)
% plot(nozzle.xPosition,nozzle.Text)
% plot(nozzle.xPosition,freestream.T*ones(length(nozzle.flow.T)))
% title('Temperature profiles')
% xlabel('Axial position (m)')
% ylabel('Temperature (K)')
% legend('T', 'T_{stag}', 'T_{w,int}', 'T_{w,ext}','T_{\infty}','Location','EastOutside');

% format plots to look nice
formatPlot;


