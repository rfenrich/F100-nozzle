% Drive the turbofanF100 function and make some interesting plots.
% SI units are used for every parameter, except for altitude (in
% feet), and sfc (in lb/lbf/hr). Sorry for mixing units.
%
% Rick Fenrich 7/23/15, modified 10/15/15

% =========================== SET INPUTS =================================

<<<<<<< HEAD
altitude = 35000; % in feet
mach = 0.65;
=======
altitude = 15000; % in feet
mach = 0.6;

% =================== SET ERROR TOLERANCE RANGES =========================
% Set error tolerances for various iterations and solvers
error.betweenIterations.inletMach = 1e-10;
error.solver.inletMach = 1e-8;
error.betweenIterations.exitTemp = 1e-6;
error.solver.apparentThroatLocation = 1e-6;
error.solver.M2relative = 1e-10;
error.solver.M2absolute = 1e-10;
error.dMdxDenominator = 4; % this is not an error tolerance, rather it is used to set the slope of dMdx in the transonic regime
>>>>>>> temporary_sensitivity_study

% ======================= INITIALIZE CONTROLS ============================
% If a control is set to zero then turbofanF100.m will
% assume a typical value for that parameter; otherwise, it will use the
% parameter value that the user provides.

% -------------------------- ENGINE CONTROLS -----------------------------

control.bypassRatio = 0;
control.f = 0;
control.fan.PstagRatio = 0;
control.fan.efficiency.polytropic = 0;
control.compressor.efficiency.polytropic = 0;
control.compressor.overallPressureRatio = 0;
control.burner.PstagRatio = 0;
control.burner.efficiency = 0;
control.turbine.efficiency.polytropic = 0;
control.turbine.efficiency.shaft = 0;

% ---------------------- NOZZLE GEOMETRY CONTROLS ------------------------

control.nozzle.geometry.shape = 'spline';
control.nozzle.geometry.length = 1;
control.nozzle.geometry.xThroat = 0.33;

if(strcmp(control.nozzle.geometry.shape,'spline'))
    % To parameterize using a spline, the following must be provided:
    % nozzle.spline.seed = either a shape already defined in the
    % nozzleGeometry.m file or an array of the form [x; y] where [x,y]
    % denote the location of the control points with the origin being at
    % the center of the inlet area
    % nozzle.spline.nControlPoints = number of control points
    % nozzle.spline.controlPointSpacing = either 'regular' where control
    % points will be evenly spaced or a vector giving the x-location
    % nozzle.spline.slopes = 1x2 array; 1st argument is slope of inlet,
    % 2nd argument is slope of outlet
    control.nozzle.geometry.spline.seed = 'linear'; %[0, 0.3255; 0.33, 0.2783; 1, 0.3293]';
    control.nozzle.geometry.spline.nControlPoints = 3;
    control.nozzle.geometry.spline.controlPointSpacing = [0 control.nozzle.geometry.xThroat control.nozzle.geometry.length]'; % 'regular';
    control.nozzle.geometry.spline.slopes = [0, 0];
end

control.nozzle.inlet.Abypass2Acore = 0;
control.nozzle.inlet.D = 0; % m
control.nozzle.throat.A = 0; % m^2

control.nozzle.geometry.Ainlet2Athroat = 1.368;
control.nozzle.geometry.Aexit2Athroat = 1.4;

% ----------------------- RUN SIMPLE TEST CASE ---------------------------
tic
[ thrust, sfc, thermalEfficiency, engine ] = turbofanF100( altitude, mach, control, error );
toc
fprintf('\n')
fprintf('thrust: %f N\n',thrust.total);
fprintf('sfc: %f lb/lbf/hr\n',sfc);
fprintf('mass flow rate: %f kg/s\n',engine.nozzle.massFlowRate);
fprintf('thermal efficiency: %f\n',thermalEfficiency);
fprintf('nozzle is: %s\n',engine.nozzle.status);
fprintf('nozzle Pstag ratio: %f\n',engine.nozzle.PstagRatio);
fprintf('nozzle Tstag ratio: %f\n',engine.nozzle.TstagRatio);
