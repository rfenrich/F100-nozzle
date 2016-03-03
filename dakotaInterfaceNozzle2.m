%% -----------------------
%%   Dakota interface script
%%   for nozzleIdeal() , nozzleNonIdeal(), and nozzleCFD()
%% -----------------------

% Rick Feinrich & Victorien Menier, last update: Feb 2016

clear all; close all;

fidelity         = 'low';
nozzle.meshSize  = 'coarse'; % 'coarse', 'medium', or 'fine'
nozzle.governing = 'euler'; % 'euler' or 'rans'

% ---------------------------
% ---  Get param file name
% ---------------------------

path = textscan(pwd,'%s','delimiter','.');
mytag = NaN;
if(numel(path{:})>1)
    mytag_temp = path{:}(end);
    mytag = str2num(mytag_temp{:});
end

if(isnan(mytag))
    parameters_file = 'params.in';
    results_file = 'results.out';
else
    parameters_file = sprintf('params.in.%d',mytag);
    results_file = sprintf('results.out.%d',mytag);
end

% ---------------------------
% --- Read dakota input file
% --- params.in (or params.in.%d )
% ---------------------------

% --- Read input parameters and update default values
params=[];
params = dakotaParseParams(parameters_file);

% ---------------------------
% --- Set default nozzle
% ---------------------------

if ( ~isnan(params.hInf) )
	nozzle.hInf = params.hInf;
	fprintf('  -- Updated input hInf = %f\n', params.hInf);
else
	nozzle.hInf = 500;
end

% --- Material properties
nozzle.wall.k                     = 8.6;    % W/m*K, thermal conductivity of nozzle wall
nozzle.wall.coeffThermalExpansion = 2.3e-6; % 1/K, coefficient of thermal expansion of nozzle wall
nozzle.wall.E                     = 80e9;   % Pa, elastic modulus of nozzle wall
nozzle.wall.poissonRatio          = 0.3;    % Poisson ratio of nozzle wall

% ---------------------------
% --- Set error tolerance ranges
% ---------------------------

error.betweenIterations.exitTemp    = 1e-8;
error.solver.apparentThroatLocation = 1e-6;
error.solver.M2relative             = 1e-10;
error.solver.M2absolute             = 1e-10;
error.dMdxDenominator               = 4; % this is not an error tolerance, rather it is used to set the slope of dMdx in the transonic regime

% ---------------------------
% --- Set nozzle geometry
% ---------------------------

nozzle.geometry.shape = 'B-spline'; % options include 'linear', 'spline', 'B-spline', and 'B-spline-mex'

% --- B-spline geometry defined by knots vector and coefficients matrix
nozzle.geometry.bSpline.knots = [0 0 0 0 1:12 13 13 13 13]';
nozzle.geometry.bSpline.coefs = [0.0000 0.0000 0.1500 0.1700 0.1900 0.2124 0.2269 0.2734 0.3218 0.3343 0.3474 0.4392 0.4828 0.5673 0.6700 0.6700;
                                   0.3255 0.3255 0.3255 0.3255 0.3255 0.3238 0.2981 0.2817 0.2787 0.2790 0.2804 0.2936 0.2978 0.3049 0.3048 0.3048];

% --- Determine nozzle throat
nozzle.geometry.bSpline.degree = length(nozzle.geometry.bSpline.knots) - length(nozzle.geometry.bSpline.coefs) - 1;
if(nozzle.geometry.bSpline.degree == 2)
    [xThroat, yThroat] = BsplineGeometry(0, 'throat', nozzle.geometry.bSpline.knots, nozzle.geometry.bSpline.coefs);
elseif(nozzle.geometry.bSpline.degree == 3)
    [xThroat, yThroat] = BsplineGeometry3(0, 'throat', nozzle.geometry.bSpline.knots, nozzle.geometry.bSpline.coefs);        
end
nozzle.geometry.xThroat = xThroat;

% --- Define other geometry parameters
nozzle.geometry.length         = nozzle.geometry.bSpline.coefs(1,end);
nozzle.inlet.D                 = nozzle.geometry.bSpline.coefs(2,1)*2;


nozzle.geometry.Ainlet2Athroat = (nozzle.inlet.D/2)^2/yThroat^2;
nozzle.geometry.Aexit2Athroat  = (nozzle.geometry.bSpline.coefs(2,end))^2/yThroat^2;

% --- Define wall geometry
nozzle.wall.shape = 'piecewise-linear';

% seed =  array of the form [x; y] where [x,y] denote the location of the control points
%         with the origin being at the center of the inlet area

if ( ~isnan(params.t1) )
	t1 = params.t1;
	fprintf('  -- Updated input t1 = %f\n', params.t1);
else
	t1 = 0.01;
end

if ( ~isnan(params.t2) )
	t2 = params.t2;
	fprintf('  -- Updated input t2 = %f\n', params.t2);
else
	t2 = 0.01;
end

if ( ~isnan(params.t3) )
	t3 = params.t3;
	fprintf('  -- Updated input t1 = %f\n', params.t3);
else
	t3 = 0.01;
end

nozzle.wall.seed = [0, t1; 
                    nozzle.geometry.xThroat, t2; 
                    nozzle.geometry.length, t3];
%  breaks = vector giving location of breaks in piecewise function
nozzle.wall.breaks = [0;
                      nozzle.geometry.xThroat;
                      nozzle.geometry.length];

% -----------------------------
% --- Set freestream conditions
% -----------------------------

fluid.gam = 1.4;    % ratio of specific heats
fluid.R   = 287.06; % J/kg-K, specific gas constant

% Flight regime 

if ( ~isnan(params.alt) )
	altitude = params.alt;
	fprintf('  -- Updated input altitude = %f\n', params.alt);
else
	altitude = 0;
end

if ( ~isnan(params.mach) )
	mach = params.mach;
	fprintf('  -- Updated input mach = %f\n', params.mach);
else
	mach = 0.01;
end

if ( ~isnan(params.Tstag) )
	nozzle.inlet.Tstag = params.Tstag;
	fprintf('  -- Updated input Tstag = %f\n', params.Tstag);
else
	nozzle.inlet.Tstag = 888.3658;
end

if ( ~isnan(params.Pstag) )
	nozzle.inlet.Pstag = params.Pstag;
	fprintf('  -- Updated input Pstag = %f\n', params.Pstag);
else
	nozzle.inlet.Pstag = 3.0550e5;
end


atm = StndAtm(altitude*0.3048,'SI'); % obtain standard atmosphere characteristics

freestream.P = atm.P;  % Pa, atmospheric pressure
freestream.T = atm.T;  % K, atmospheric temperature
freestream.M = mach;
freestream.U = freestream.M*sqrt(fluid.gam*fluid.R*freestream.T);

% ---------------------------
% --- Run nozzle computation
% ---------------------------

if(strcmp(fidelity,'low'))
	[ nozzle ] = nozzleNonIdeal( fluid, freestream, nozzle, error );
elseif(strcmp(fidelity,'medium'))
	[ nozzle ] = nozzleCFD( fluid, freestream, nozzle, error );
else
	error('Incorrect keyword for specifying the fidelity level.');
end

% ---------------------------
% --- Output functional of interest : thrust
% ---------------------------

fprintf('  -- Thrust     = %f\n', nozzle.netThrust);


fid = fopen(results_file,'w');
fprintf(fid,'%.16e thrust\n',nozzle.netThrust);
fprintf(fid,'%.16e maxStress\n',max(nozzle.hoopStress));
fprintf(fid,'%.16e maxTemp\n',nozzle.exit.T);
fprintf(fid,'%.16e volume\n',nozzle.geometry.volume);
fclose(fid);
