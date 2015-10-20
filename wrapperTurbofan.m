function [thrust, sfc, thermalEfficiency] = wrapperTurbofan(x)

% Normalize the parameters
% 0.57 <= control.bypassRatio <= 0.63 (nominal 0.6)
% 2.91 <= control.fan.PstagRatio <= 3.21 (nominal 3.06)
% 0.82 <= control.fan.efficiency.polytropic <= 0.86 
% 24 <= control.compressor.overallPressureRatio <= 25 (nominal 24.5)
% 0.84 <= control.compressor.efficiency.polytropic <= 0.9
% 0.92 <= control.burner.PstagRatio <= 0.98
% 0.94 <= control.burner.efficiency <= 0.99
% 0.83 <= control.turbine.efficiency.polytropic <= 0.89
% 0.95 <= control.turbine.efficiency.shaft <= 0.99
% 0.15 <= control.nozzle.inlet.Abypass2Acore <= 0.4 (nominal 0.296)

lb = [0.57 2.91 0.82 24 0.84 0.92 0.94 0.83 0.95 0.15];
ub = [0.63 3.21 0.86 25 0.9 0.98 0.99 0.89 0.99 0.4];
x0 = 0.5*(x+1).*(ub-lb) + lb;

% Set parameters of turbofan script
altitude = 35000; % in feet
mach = 0.9;

control.bypassRatio = x0(1);
control.f = 0;
control.fan.PstagRatio = x0(2);
control.fan.efficiency.polytropic = x0(3);
control.compressor.efficiency.polytropic = x0(5);
control.compressor.overallPressureRatio = x0(4);
control.burner.PstagRatio = x0(6);
control.burner.efficiency = x0(7);
control.turbine.efficiency.polytropic = x0(8);
control.turbine.efficiency.shaft = x0(9);
control.nozzle.inlet.Abypass2Acore = x0(10);
control.nozzle.inlet.D = 0; % m
control.nozzle.throat.A = 0; % m^2

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

control.nozzle.inlet.D = 0; % m
control.nozzle.throat.A = 0; % m^2

control.nozzle.geometry.Ainlet2Athroat = 1.368;
control.nozzle.geometry.Aexit2Athroat = 1.4;

% Set error tolerances for various iterations and solvers
error.betweenIterations.inletMach = 1e-10;
error.solver.inletMach = 1e-8;
error.betweenIterations.exitTemp = 1e-6;
error.solver.apparentThroatLocation = 1e-6;
error.solver.M2relative = 1e-10;
error.solver.M2absolute = 1e-10;
error.dMdxDenominator = 4; % this is not an error tolerance, rather it is used to set the slope of dMdx in the transonic regime

[thrust, sfc, thermalEfficiency, ~] = turbofanF100( altitude, mach, control, error );


