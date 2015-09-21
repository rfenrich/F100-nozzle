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
control.nozzle.Ainlet2Athroat = 0;
control.nozzle.Aexit2Athroat = 0;

control.nozzle.Ainlet2Athroat = 1.368;
control.nozzle.Aexit2Athroat = 1.4;

[thrust, sfc, thermalEfficiency, ~] = turbofanF100( altitude, mach, control );


