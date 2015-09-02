% Drive the turbofanF100 function and make some interesting plots.
% SI units are used for every parameter, except for altitude (in
% feet), and sfc (in lb/lbf/hr). Sorry for mixing units.
%
% Rick Fenrich 7/23/15, modified 7/31/15

% =========================== SET INPUTS =================================

altitude = 35000; % in feet
mach = 0.9;

% ======================= INITIALIZE CONTROLS ============================
% If a control is set to zero then turbofanF100.m will
% assume a typical value for that parameter; otherwise, it will use the
% parameter value that the user provides.

control.bypassRatio = 0;
control.f = 0;
control.fan.PstagRatio = 0;
control.fan.efficiency.polytropic = 0;
control.fan.exit.M = 0;
control.compressor.efficiency.polytropic = 0;
control.compressor.overallPressureRatio = 0;
control.burner.PstagRatio = 0;
control.burner.efficiency = 0;
control.turbine.efficiency.polytropic = 0;
control.turbine.efficiency.shaft = 0;
control.nozzle.inlet.Abypass2Acore = 0;
control.nozzle.inlet.D = 0; % m
control.nozzle.throat.A = 0; % m^2
control.nozzle.Ainlet2Athroat = 0;
control.nozzle.Aexit2Athroat = 0;

control.nozzle.Ainlet2Athroat = 1.368;
control.nozzle.Aexit2Athroat = 1.4;
control.fan.exit.M = 0.8;

% ----------------------- RUN SIMPLE TEST CASE ---------------------------
tic
[ thrust, sfc, thermalEfficiency, engine ] = turbofanF100( altitude, mach, control );
toc
fprintf('\n')
fprintf('thrust: %f N\n',thrust.total);
fprintf('sfc: %f lb/lbf/hr\n',sfc);
fprintf('mass flow rate: %f kg/s\n',engine.nozzle.massFlowRate);
fprintf('thermal efficiency: %f\n',thermalEfficiency);
fprintf('nozzle is: %s\n',engine.nozzle.status);