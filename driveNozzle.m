% ========= Driver script for nozzleIdeal.m and nozzleNonIdeal.m =========
%
% Provides inputs to the nozzleIdeal and nozzleNonIdeal functions, runs
% them, and produces plots comparing ideal and non-ideal nozzle behavior,
% as well as non-ideal nozzle heating and friction characteristics.
%
% Rick Fenrich 7/29/15 modified 2/17/16

% ========================== INPUT PARAMETERS ============================

mission = 3; % mission number for which certain input parameters are defined below
fidelity = 'low'; % 'low' (quasi-1D flow), 'med' (Euler), or 'high' (TBD)

% ============= Options for mid-fidelity (used only if fidelity=='med')

nozzle.meshSize  = 'coarse'; % 'coarse', 'medium', or 'fine'
nozzle.governing = 'euler'; % 'euler' or 'rans'

% --------------------- SET HEAT TRANSFER PARAMS -------------------------
nozzle.hInf = 500; % W/m^2-K, heat transfer coefficient from external nozzle wall to environment

% --------------------- SET MATERIAL PROPERTIES --------------------------
% SEPCARBINOX A500, a ceramic matrix composite
nozzle.wall.k = 8.6; % W/m*K, thermal conductivity of nozzle wall
nozzle.wall.coeffThermalExpansion = 2.3e-6; % 1/K, coefficient of thermal expansion of nozzle wall
nozzle.wall.E = 80e9; % Pa, elastic modulus of nozzle wall
nozzle.wall.poissonRatio = 0.3; % Poisson ratio of nozzle wall

% ------------------------ SET FLIGHT REGIME -----------------------------
% Define input parameters that will change based on flight regime:
if(mission == 1) % static sea-level thrust case
    altitude = 0;
    mach = 0.01;
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
    mach = 0.01;
    nozzle.inlet.Tstag = 900;
    nozzle.inlet.Pstag = 1.3e5;
elseif(mission == 5) % subsonic flow in nozzle
    altitude = 0;
    mach = 0.01;
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
nozzle.geometry.shape = 'B-spline-mex'; % options include 'linear', 'spline', 'B-spline', and 'B-spline-mex'
if(strcmp(nozzle.geometry.shape,'spline')) 
    % To parameterize using a cubic spline, the following must be provided:
    % spline.seed = either a shape already defined in the nozzleGeometry.m 
    %   file (such as 'linear') or an array of the form [x,y] where [x,y] 
    %   denote the location of the control points with the origin being at 
    %   the center of the inlet area
    % spline.slopes = 1x2 array; 1st argument is slope of inlet, 2nd 
    %   argument is slope of outlet
    % spline.breaks = a vector giving the x-location of control points 
    %   (only necessary if spline.seed is NOT a matrix of coordinates)

    nozzle.geometry.spline.seed = ...%'linear';
                                [0, 0.3255; 
                                 0.33, 0.2783; 
                                 1, 0.3293];
    nozzle.geometry.spline.slopes = [0, 0];
       
    % Placeholder for throat coordinates (will be recalculated later)
    nozzle.geometry.xThroat = 0.33; % m, location of throat from inlet
    nozzle.geometry.yThroat = 0.2783;
    
    % Spline breaks only need to be specified if coordinates are not given
    % in the nozzle.geometry.spline.seed
    nozzle.geometry.spline.breaks = [0; nozzle.geometry.xThroat; nozzle.geometry.length];

    % Recalculate area ratios/inlet geometry if necessary from the spline seed:
    nozzle.geometry.length = nozzle.geometry.spline.seed(end,1); % m
    nozzle.inlet.D = (nozzle.geometry.spline.seed(1,2))*2;
    nozzle.geometry.Ainlet2Athroat = (nozzle.inlet.D/2)^2/(nozzle.geometry.yThroat)^2; % area ratio of inlet to throat
    nozzle.geometry.Aexit2Athroat = (nozzle.geometry.spline.seed(end,2))^2/(nozzle.geometry.yThroat)^2; % area ratio of exit to throat
    
elseif(strcmp(nozzle.geometry.shape,'B-spline') || strcmp(nozzle.geometry.shape,'B-spline-mex'))
    % B-spline geometry defined by knots vector and coefficients matrix
	nozzle.geometry.bSpline.knots = [0 0 0 0 1:12 13 13 13 13]';
	nozzle.geometry.bSpline.coefs = [0.0000 0.0000 0.1500 0.1700 0.1900 0.2124 0.2269 0.2734 0.3218 0.3343 0.3474 0.4392 0.4828 0.5673 0.6700 0.6700;
                                     0.3255 0.3255 0.3255 0.3255 0.3255 0.3238 0.2981 0.2817 0.2787 0.2790 0.2804 0.2936 0.2978 0.3049 0.3048 0.3048];

    % Determine nozzle throat
    nozzle.geometry.bSpline.degree = length(nozzle.geometry.bSpline.knots) - length(nozzle.geometry.bSpline.coefs) - 1;
    if(nozzle.geometry.bSpline.degree == 2)
        [xThroat, yThroat] = BsplineGeometry(0, 'throat', nozzle.geometry.bSpline.knots, nozzle.geometry.bSpline.coefs);
    elseif(nozzle.geometry.bSpline.degree == 3)
        [xThroat, yThroat] = BsplineGeometry3(0, 'throat', nozzle.geometry.bSpline.knots, nozzle.geometry.bSpline.coefs);        
    end
    nozzle.geometry.xThroat = xThroat;

    % Define other geometry parameters
    nozzle.geometry.length = nozzle.geometry.bSpline.coefs(1,end);
    nozzle.inlet.D = nozzle.geometry.bSpline.coefs(2,1)*2;
    nozzle.geometry.Ainlet2Athroat = (nozzle.inlet.D/2)^2/yThroat^2;
    nozzle.geometry.Aexit2Athroat = (nozzle.geometry.bSpline.coefs(2,end))^2/yThroat^2;
    
end

% -------------------- SET NOZZLE WALL GEOMETRY --------------------------
nozzle.wall.shape = 'piecewise-linear';
% nozzle.wall.seed = an array of the form [x; y] where [x,y] denote the 
%   location of the control points with the origin being at the center of 
%   the inlet area
% nozzle.wall.breaks = vector giving location of breaks in piecewise
%   function
nozzle.wall.seed = [0, 0.01; 
                    nozzle.geometry.xThroat, 0.01; 
                    nozzle.geometry.length, 0.01];
nozzle.wall.breaks = [0;
                      nozzle.geometry.xThroat;
                      nozzle.geometry.length];

% --------------------- FLUID INPUT PROPERTIES ---------------------------
fluid.gam = 1.4; % ratio of specific heats
fluid.R = 287.06; % J/kg-K, specific gas constant

% ====================== RUN NOZZLE CALCULATIONS =========================

% ----------------- CALCULATE FREESTREAM PROPERTIES ----------------------
atm = StndAtm(altitude*0.3048,'SI'); % obtain standard atmosphere characteristics
freestream.P = atm.P; % Pa, atmospheric pressure
freestream.T = atm.T; % K, atmospheric temperature
freestream.M = mach;
freestream.U = freestream.M*sqrt(fluid.gam*fluid.R*freestream.T);

% ----------------- PRINT SOME USEFUL DATA TO SCREEN ---------------------
fprintf('Mission %i: %i ft altitude at %0.2f Mach\n',mission,altitude,mach);
if(strcmp(nozzle.geometry.shape,'B-spline') || strcmp(nozzle.geometry.shape,'B-spline-mex'))
    fprintf('Geometry: %s (degree %i)\n',nozzle.geometry.shape,nozzle.geometry.bSpline.degree);
else
    fprintf('Geometry: %s \n',nozzle.geometry.shape);
end
fprintf('Wall Geo: %s\n',nozzle.wall.shape);

% ------------------------- RUN CALCULATIONS -----------------------------
if(strcmp(fidelity,'low'))
	[ nozzleI ] = nozzleIdeal( fluid, freestream, nozzle, error );
    tic;
	[ nozzle ] = nozzleNonIdeal( fluid, freestream, nozzle, error );
    TimeToExec = toc;
elseif(strcmp(fidelity,'med'))
    tic;
    [ nozzle ] = nozzleCFD( fluid, freestream, nozzle, error);
    TimeToExec = toc;
else
    error('Incorrect fidelity specified.');
end

% ============================ OUTPUT DATA ===============================

fprintf('====== RESULTS ======\n');
fprintf('Time to execute: %0.3f\n',TimeToExec);
fprintf('Pstag ratio: %0.4f\n',nozzle.exit.Pstag/nozzle.inlet.Pstag);
fprintf('Tstag ratio: %0.4f\n',nozzle.exit.Tstag/nozzle.inlet.Tstag);

% Output other data useful for CFD simulation:
fprintf('Freestream Mach: %0.2f\n',mach);
fprintf('Freestream Static Pressure: %0.4f Pa\n',freestream.P);
fprintf('Freestream Static Temperature: %0.4f K\n',freestream.T);
fprintf('Inlet Stagnation Pressure: %0.4f Pa\n',nozzle.inlet.Pstag);
fprintf('Inlet Stagnation Temperature: %0.4f K\n',nozzle.inlet.Tstag);
fprintf('Outlet Static Pressure: %0.4f Pa\n',nozzle.exit.P);

fprintf('Nozzle Volume: %0.4f cm^3\n',nozzle.geometry.volume*100^3);
fprintf('Nozzle Est. Thrust: %0.4f N\n',nozzle.netThrust);

% ============================== PLOTTING ================================
% Comment or uncomment the plotting below to make plots

% Manage colors and linestyles for plots:
set(groot,'defaultAxesColorOrder',[0 0 0; 1 0 0; 1 0.25 0; 0 1 0; 0 0 1])
%set(groot,'defaultAxesLineStyleOrder','-|--|:');

% % -------------------- PLOT VARIOUS DATA IN 1 FIGURE ---------------------
% figure
% subplot(2,3,1); hold on
% plot(nozzle.xPosition,nozzle.geometry.D/2)
% plot(nozzle.xPosition,nozzle.geometry.D/2+nozzle.wall.t);
% plot(nozzle.xPosition,-nozzle.geometry.D/2);
% plot(nozzle.xPosition,-nozzle.geometry.D/2-nozzle.wall.t);
% title('Geometry')
% axis equal
% 
% subplot(2,3,2); hold on
% plot(nozzleI.xPosition,nozzleI.flow.Re)
% plot(nozzle.xPosition,nozzle.flow.Re)
% title('Re')
% legend('ideal');
% 
% subplot(2,3,3); hold on
% plot(nozzleI.xPosition,nozzleI.flow.M)
% plot(nozzle.xPosition,nozzle.flow.M)
% title('Mach Number')
% legend('ideal');
% 
% subplot(2,3,4); hold on
% plot(nozzleI.xPosition,nozzleI.flow.P)
% plot(nozzle.xPosition,nozzle.flow.P)
% plot(nozzleI.xPosition,nozzleI.flow.Pstag)
% plot(nozzle.xPosition,nozzle.flow.Pstag)
% plot(nozzle.xPosition,freestream.P*ones(length(nozzle.xPosition),1))
% title('Pressure (Pa)')
% legend('Static (ideal)','Static','Stag (ideal)','Stag','\infty');
% 
% subplot(2,3,5); hold on
% plot(nozzleI.xPosition,nozzleI.flow.T)
% plot(nozzle.xPosition,nozzle.flow.T)
% plot(nozzleI.xPosition,nozzleI.flow.Tstag)
% plot(nozzle.xPosition,nozzle.flow.Tstag)
% plot(nozzle.xPosition,freestream.T*ones(length(nozzle.xPosition),1),'k-')
% title('Temperature (K)')
% legend('Static (ideal)','Static','Stag (ideal)','Stag','\infty');
% 
% subplot(2,3,6); hold on
% plot(nozzleI.xPosition,nozzleI.flow.density)
% plot(nozzle.xPosition,nozzle.flow.density)
% title('Density (kg/m^3)')
% legend('ideal');
% 
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
% % -------------------------- PLOT MAX STRESS -----------------------------
% figure
% plot(nozzle.xPosition,nozzle.maxStress/1e6)
% xlabel('Axial position (m)')
% title('Principal Stress (MPa)')
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
% 
% % format plots to look nice
% formatPlot;

