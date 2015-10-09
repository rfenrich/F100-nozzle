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

altitude = NaN;
mach = NaN;
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
    elseif(strcmp(label_i,'bypass'))
        control.bypassRatio = var_i;
    elseif(strcmp(label_i,'fanPstag'))
        control.fan.PstagRatio = var_i;
    elseif(strcmp(label_i,'fanEff'))
        control.fan.efficiency.polytropic = var_i;
    elseif(strcmp(label_i,'compressEff'))
        control.compressor.efficiency.polytropic = var_i;
    elseif(strcmp(label_i,'compressPratio'))
        control.compressor.overallPressureRatio = var_i;
    elseif(strcmp(label_i,'burnerPstag'))
        control.burner.PstagRatio = var_i;
    elseif(strcmp(label_i,'burnerEff'))
        control.burner.efficiency = var_i;
    elseif(strcmp(label_i,'turbineEffPoly'))
        control.turbine.efficiency.polytropic = var_i;
    elseif(strcmp(label_i,'turbineEffShaft'))
        control.turbine.efficiency.shaft = var_i;
    elseif(strcmp(label_i,'Abypass2Acore'))
        control.nozzle.inlet.Abypass2Acore = var_i;
    elseif(strcmp(label_i,'nozzleInletD'))
        control.nozzle.inlet.D = var_i;
    elseif(strcmp(label_i,'Ainlet2Athroat'))
        control.nozzle.Ainlet2Athroat = var_i;
    elseif(strcmp(label_i,'Aexit2Athroat'))
        control.nozzle.Aexit2Athroat = var_i;
    else
        disp('Unknown variable!');
    end
end
fclose(fid);

% ======================= INITIALIZE CONTROLS ============================
% If a control is set to zero then turbofanF100.m will
% assume a typical value for that parameter; otherwise, it will use the
% parameter value that the user provides.

if(isnan(altitude))
    altitude = 0;
end

if(isnan(mach))
    mach = 0;
end

% ----------------------- RUN turbofan simulation ---------------------------
[ thrust, sfc, thermalEfficiency, engine ] = turbofanF100( altitude, mach, control );

% Write results to file for Dakota:
% fid = fopen(results_file,'w');
% fprintf(fid,'%f sfc',sfc);
% fclose(fid);


fid = fopen(results_file,'w');
fprintf(fid,'%f thrust\n',thrust.total);
fprintf(fid,'%f sfc\n',sfc);
fprintf(fid,'%f massFlowRate\n',engine.nozzle.massFlowRate);
fprintf(fid,'%f thermalEfficiency\n',thermalEfficiency);
fclose(fid);

% Turn singular matrix warnings back on.
warning('on','MATLAB:nearlySingularMatrix');