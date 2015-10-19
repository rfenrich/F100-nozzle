% Drive the turbofanF100 function and determine sensitivity of common
% outputs (thrust, sfc, eta, etc.) to each of the required error tolerances
% in the turbofanF100 and nonIdealNozzle functions.
%
% Rick Fenrich 10/7/15

% =========================== SET INPUTS =================================

altitude = 0; % in feet
mach = 0;

% =================== SET ERROR TOLERANCE RANGES =========================
inletMachIterationError = [1e-6 1e-8 1e-10 1e-12];
inletMachSolverError = [1e-4 1e-6 1e-8 1e-10];
exitTempIterationError = [1e-2 1e-4 1e-6 1e-8 1e-10 1e-12];
apparentThroatSolverError = [1e-2 1e-4 1e-6 1e-8];
M2SolverRelError = [1e-6 1e-8 1e-10 1e-12 1e-14];
M2SolverAbsError = [1e-6 1e-8 1e-10 1e-12 1e-14];
dMdxDenominator = [3 4 5];

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

% Open file for writing data to
fid = fopen('error_sensitivity.txt','a');
fprintf(fid,'Error sensitivity for F100-PW-220 model\n');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'altitude: %f, mach: %f\n',altitude,mach);

% Run nominal case
fprintf(fid,'Nominal case\n');
nominalIndices = [3 3 4 3 3 3 2];
fprintf(fid,'inlet Mach iter error: %e\n',inletMachIterationError(nominalIndices(1)));
fprintf(fid,'inlet Mach solver error: %e\n',inletMachSolverError(nominalIndices(2)));
fprintf(fid,'exit T iter error: %e\n',exitTempIterationError(nominalIndices(3)));
fprintf(fid,'apparent throat location solver error: %e\n',apparentThroatSolverError(nominalIndices(4)));
fprintf(fid,'M2 rel & abs solver error: %e\n',M2SolverRelError(nominalIndices(5)));
fprintf(fid,'dMdx denominator: %f\n',dMdxDenominator(nominalIndices(7)));

error.betweenIterations.inletMach = inletMachIterationError(nominalIndices(1));
error.solver.inletMach = inletMachSolverError(nominalIndices(2));
error.betweenIterations.exitTemp = exitTempIterationError(nominalIndices(3));
error.solver.apparentThroatLocation = apparentThroatSolverError(nominalIndices(4));
error.solver.M2relative = M2SolverRelError(nominalIndices(5));
error.solver.M2absolute = M2SolverAbsError(nominalIndices(6));
error.dMdxDenominator = dMdxDenominator(nominalIndices(7));
nominalError = error;

tic;
[ thrust, sfc, thermalEfficiency, engine ] = turbofanF100( altitude, mach, control, error );
time = toc;

nominal.thrust = thrust;
nominal.sfc = sfc;
nominal.thermalEfficiency = thermalEfficiency;
nominal.massFlowRate = engine.nozzle.massFlowRate;
nominal.PstagRatio = engine.nozzle.PstagRatio;
nominal.TstagRatio = engine.nozzle.TstagRatio;
nominal.time = time;

fprintf(fid,'error \t\t time \t thrust \t sfc \t mdot \t eta \t PstagRatio \t TstagRatio \n');
fprintf(fid,'%s \t %e \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \n','nominal', ...
    0, time, thrust.total, sfc, engine.nozzle.massFlowRate, thermalEfficiency, ...
    engine.nozzle.PstagRatio, engine.nozzle.TstagRatio);

% Run all non-nominal cases

for ii = 1:length(inletMachIterationError)
    error = nominalError;
    error.betweenIterations.inletMach = inletMachIterationError(ii);
    tic;
    [ thrust, sfc, thermalEfficiency, engine ] = turbofanF100( altitude, mach, control, error );
    time = toc;
    fprintf(fid,'%s \t %e \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \n','Mi iter', ...
    inletMachIterationError(ii), time, thrust.total, sfc, ...
    engine.nozzle.massFlowRate, thermalEfficiency, ...
    engine.nozzle.PstagRatio, engine.nozzle.TstagRatio);
end

for ii = 1:length(inletMachSolverError)
    error = nominalError;
    error.solver.inletMach = inletMachSolverError(ii);
    tic;
    [ thrust, sfc, thermalEfficiency, engine ] = turbofanF100( altitude, mach, control, error );
    time = toc;    
    fprintf(fid,'%s \t %e \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \n','Mi solv', ...
    inletMachSolverError(ii), time, thrust.total, sfc, ...
    engine.nozzle.massFlowRate, thermalEfficiency, ...
    engine.nozzle.PstagRatio, engine.nozzle.TstagRatio);
end

for ii = 1:length(exitTempIterationError)
    error = nominalError;
    error.betweenIterations.exitTemp = exitTempIterationError(ii);
    tic;
    [ thrust, sfc, thermalEfficiency, engine ] = turbofanF100( altitude, mach, control, error );
    time = toc;    
    fprintf(fid,'%s \t %e \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \n','Te iter', ...
    exitTempIterationError(ii), time, thrust.total, sfc, ...
    engine.nozzle.massFlowRate, thermalEfficiency, ...
    engine.nozzle.PstagRatio, engine.nozzle.TstagRatio);
end

for ii = 1:length(apparentThroatSolverError)
    error = nominalError;
    error.solver.apparentThroatLocation = apparentThroatSolverError(ii);
    tic;
    [ thrust, sfc, thermalEfficiency, engine ] = turbofanF100( altitude, mach, control, error );
    time = toc;    
    fprintf(fid,'%s \t %e \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \n','th solv', ...
    apparentThroatSolverError(ii), time, thrust.total, sfc, ...
    engine.nozzle.massFlowRate, thermalEfficiency, ...
    engine.nozzle.PstagRatio, engine.nozzle.TstagRatio);
end

for ii = 1:length(M2SolverRelError)
    error = nominalError;
    error.solver.M2relative = M2SolverRelError(ii);
    error.solver.M2absolute = M2SolverAbsError(ii);
    tic;
    [ thrust, sfc, thermalEfficiency, engine ] = turbofanF100( altitude, mach, control, error );
    time = toc;    
    fprintf(fid,'%s \t %e \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \n','M2 solv', ...
    M2SolverRelError(ii), time, thrust.total, sfc, ...
    engine.nozzle.massFlowRate, thermalEfficiency, ...
    engine.nozzle.PstagRatio, engine.nozzle.TstagRatio);
end

for ii = 1:length(dMdxDenominator)
    error = nominalError;
    error.dMdxDenominator = dMdxDenominator(ii);
    tic;
    [ thrust, sfc, thermalEfficiency, engine ] = turbofanF100( altitude, mach, control, error );
    time = toc;    
    fprintf(fid,'%s \t %e \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \t %1.12f \n','dMdxden', ...
    dMdxDenominator(ii), time, thrust.total, sfc, ...
    engine.nozzle.massFlowRate, thermalEfficiency, ...
    engine.nozzle.PstagRatio, engine.nozzle.TstagRatio);
end

fclose(fid);
