% Interface between turbofanF100 function and Dakota.
% Based off of driveTurbofanModel.m
%
% Rick Fenrich 7/23/15, modified 7/31/15
% Modified by Jason Monschke 10/9/15

% =========================== SET INPUTS =================================

clc; clear all; close all;
% Turn off singular matrix warnings from trustsolve.m (In this case they
% are expected).
warning('off','MATLAB:nearlySingularMatrix');

% Read parameters from Dakota generated parameter file:
parameters_file = 'params.in';
results_file = 'results.out';

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
    else
        disp('Unknown variable!');
    end
end
fclose(fid);

% ======================= INITIALIZE CONTROLS ============================
% If a control is set to zero then turbofanF100.m will
% assume a typical value for that parameter; otherwise, it will use the
% parameter value that the user provides.

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
control.nozzle.inlet.Abypass2Acore = 0;
control.nozzle.inlet.D = 0; % m
control.nozzle.throat.A = 0; % m^2
control.nozzle.Ainlet2Athroat = 0;
control.nozzle.Aexit2Athroat = 0;

control.nozzle.Ainlet2Athroat = 1.368;
control.nozzle.Aexit2Athroat = 1.4;

% ----------------------- RUN turbofan simulation ---------------------------
[ thrust, sfc, thermalEfficiency, engine ] = turbofanF100( altitude, mach, control );

% Write results to file for Dakota:
fid = fopen(results_file,'w');
fprintf(fid,'%f sfc',sfc);
fclose(fid);

% Turn singular matrix warnings back on.
warning('on','MATLAB:nearlySingularMatrix');