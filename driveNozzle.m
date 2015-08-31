% ========= Driver script for nozzleIdeal.m and nozzleNonIdeal.m =========
%
% Provides inputs to the nozzleIdeal and nozzleNonIdeal functions, runs
% them, and produces plots comparing ideal and non-ideal nozzle behavior,
% as well as non-ideal nozzle heating and friction characteristics.
%
% Rick Fenrich 7/29/15

% ========================== INPUT PARAMETERS ============================

mission = 4; % select the mission which defines input parameters
hInf = 500; % W/m^2-K, heat transfer coefficient from external nozzle wall to environment

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
    nozzle.inlet.Pstag = 1.4492e5;
elseif(mission == 4) % case with shock in nozzle
    altitude = 0;
    mach = 0.0;
    nozzle.inlet.Tstag = 900;
    nozzle.inlet.Pstag = 1.1e5;
end

% Other necessary input parameters
fluid.gam = 1.4;
fluid.R = 287.06;
nozzle.inlet.D = 0.651;
nozzle.Ainlet2Athroat = 1.368;
nozzle.Aexit2Athroat = 1.4;
nozzle.length = 1;
nozzle.shape = 'linear';
nozzle.xThroat = 0.33;
nozzle.xExit = nozzle.length;

atm = StndAtm(altitude*0.3048,'SI');
freestream.P = atm.P; % Pa, atmospheric pressure
freestream.T = atm.T; % K, atmospheric temperature

% ====================== RUN NOZZLE CALCULATIONS =========================

[ flow.ideal, nozzle.ideal, xPositionIdeal ] = nozzleIdeal( fluid, nozzle.inlet, freestream, nozzle);

[ flow.nonideal, nozzle.nonideal, xPosition ] = nozzleNonIdeal( fluid, nozzle.inlet, freestream, nozzle, hInf);

% ============================ OUTPUT DATA ===============================

fprintf('Pstag ratio: %f\n',flow.nonideal.Pstag(end)/flow.nonideal.Pstag(1));
fprintf('Tstag ratio: %f\n',flow.nonideal.Tstag(end)/flow.nonideal.Tstag(1));

% ============================== PLOTTING ================================
% Comment or uncomment the plotting below to make plots

% Manage colors and linestyles for plots:
set(groot,'defaultAxesColorOrder',[0 0 0; 1 0 0; 1 0.25 0; 0 1 0; 0 0 1])
%set(groot,'defaultAxesLineStyleOrder','-|--|:');

% Useful function which defines diameter of nozzle; for plotting
D = @(x) nozzleGeometry(x,nozzle.inlet.D,nozzle.length,nozzle.xThroat,nozzle.Ainlet2Athroat,nozzle.Aexit2Athroat,nozzle.shape,'D');

% -------------------- PLOT VARIOUS DATA IN 1 FIGURE ---------------------
figure
subplot(2,3,1); hold on
plot(xPosition,D(xPosition)/2)
plot(xPosition,-D(xPosition)/2)
title('Geometry')
axis equal

subplot(2,3,2); hold on
plot(xPositionIdeal,flow.ideal.Re)
plot(xPosition,flow.nonideal.Re)
title('Re')
legend2 = legend('ideal');

subplot(2,3,3); hold on
plot(xPositionIdeal,flow.ideal.M)
plot(xPosition,flow.nonideal.M)
title('Mach Number')
legend3 = legend('ideal');

subplot(2,3,4); hold on
plot(xPositionIdeal,flow.ideal.P)
plot(xPosition,flow.nonideal.P)
plot(xPositionIdeal,flow.ideal.Pstag)
plot(xPosition,flow.nonideal.Pstag)
plot(xPosition,freestream.P*ones(length(xPosition),1))
title('Pressure (Pa)')
legend4 = legend('Static (ideal)','Static','Stag (ideal)','Stag','\infty');

subplot(2,3,5); hold on
plot(xPositionIdeal,flow.ideal.T)
plot(xPosition,flow.nonideal.T)
plot(xPositionIdeal,flow.ideal.Tstag)
plot(xPosition,flow.nonideal.Tstag)
plot(xPosition,freestream.T*ones(length(xPosition),1),'k-')
title('Temperature (K)')
legend5 = legend('Static (ideal)','Static','Stag (ideal)','Stag','\infty');

subplot(2,3,6); hold on
plot(xPositionIdeal,flow.ideal.density)
plot(xPosition,flow.nonideal.density)
title('Density (kg/m^3)')
legend6 = legend('ideal');

% ------------------------ PLOT NOZZLE GEOMETRY --------------------------
figure; hold on
plot(xPosition,D(xPosition)/2)
plot(xPosition,-D(xPosition)/2)
xlabel('Axial position (m)')
title('Geometry')
axis equal

% ------------------------ PLOT REYNOLDS NUMBER --------------------------
figure; hold on
plot(xPositionIdeal,flow.ideal.Re);
plot(xPosition,flow.nonideal.Re)
xlabel('Axial position (m)')
title('Re')
legend2 = legend('ideal');

% ------------------------- PLOT MACH GEOMETRY ---------------------------
figure; hold on
plot(xPositionIdeal,flow.ideal.M)
plot(xPosition,flow.nonideal.M)
xlabel('Axial position (m)')
title('Mach Number')
legend3 = legend('ideal');

% --------------------------- PLOT PRESSURES -----------------------------
figure; hold on
plot(xPositionIdeal,flow.ideal.P)
plot(xPosition,flow.nonideal.P)
plot(xPositionIdeal,flow.ideal.Pstag)
plot(xPosition,flow.nonideal.Pstag)
plot(xPosition,freestream.P*ones(length(xPosition),1))
xlabel('Axial position (m)')
title('Pressure (Pa)')
legend4 = legend('Static (ideal)','Static','Stag (ideal)','Stag','\infty','Location','EastOutside');

% ------------------------- PLOT TEMPERATURES ----------------------------
figure; hold on
hold on
plot(xPositionIdeal,flow.ideal.T)
plot(xPosition,flow.nonideal.T)
plot(xPositionIdeal,flow.ideal.Tstag)
plot(xPosition,flow.nonideal.Tstag)
plot(xPosition,freestream.T*ones(length(xPosition),1))
xlabel('Axial position (m)','FontName','CMU Serif','FontSize',14)
title('Temperature (K)','FontName','CMU Serif','FontSize',16)
set(gca,'FontName','CMU Serif','FontSize',14)
legend5 = legend('Static (ideal)','Static','Stag (ideal)','Stag','\infty','Location','EastOutside');
set(legend5,'FontName','CMU Serif','FontSize',12)

% ---------------------------- PLOT DENSITY ------------------------------
figure; hold on
plot(xPositionIdeal,flow.ideal.density)
plot(xPosition,flow.nonideal.density)
xlabel('Axial position (m)')
title('Density (kg/m^3)')
legend6 = legend('ideal');

% --------------------------- PLOT h_f & C_f -----------------------------
figure
subplot(1,2,1); hold on
plot(xPosition,flow.nonideal.hf)
xlabel('Axial position (m)')
title('Wall convection coefficient h_f (W/m^2-K)')
subplot(1,2,2); hold on
plot(xPosition,flow.nonideal.Cf)
xlabel('Axial position (m)')
title('Friction coefficient C_f')

% --------------------- PLOT TEMPERATURE PROFILES ------------------------
figure
hold on
plot(xPosition,flow.nonideal.T)
plot(xPosition,flow.nonideal.Tstag)
plot(xPosition,nozzle.nonideal.Tw)
plot(xPosition,nozzle.nonideal.Text)
plot(xPosition,freestream.T*ones(length(flow.nonideal.T)))
title('Temperature profiles')
xlabel('Axial position (m)')
ylabel('Temperature (K)')
legend7 = legend('T', 'T_{stag}', 'T_{w,int}', 'T_{w,ext}','T_{\infty}','Location','EastOutside');

formatPlot;

