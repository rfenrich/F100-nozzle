function [varargout] = exampleWrapper(x,knots,coefs,fidelity,output,F100dir)
% Example wrapper for an optimization using nonIdealNozzle.m. This wrapper
% in particular corresponds to exampleOptDriver.m. The inputs knots, coefs,
% and objective can be defined inside this wrapper, so that exampleWrapper
% is only a function of x, the design variables.
%
% exampleOptDriver.m uses the x and y coordinates of the control points of
% a B-spline parameterized nozzle geometry as design variables to optimize
% deterministically for min weight subject to a thrust constraint.
%
% INPUTS:
% x = input vector of design variables
% knots = knots vector required for B-spline parameterization
% coefs = coefficient matrix required for B-spline parameterization
% fidelity = 'low' or 'med'; determines which nozzle calculation to run
% output = string containing name of output variable
% F100dir = string containing absolute directory of F100-nozzle code;
% useful if optimization is run in a different directory
%
% OUTPUTS:
% varargout = output of wrapper specified by objective string in inputs
%
% Rick Fenrich 2/17/16 updated 5/17/16

% ====================== ASSIGN DESIGN VARIABLES =========================

% Redefine B-spline geometry using design variables
coefs(1,6:9) = x(1:4);
coefs(1,10) = x(4);
coefs(1,11:17) = x(5:11);
coefs(1,18) = x(11);
coefs(2,6:9) = x(12:15);
coefs(2,10) = x(15);
coefs(2,11) = x(15);
coefs(2,12:17) = x(16:21);
coefs(2,18) = x(21);
nozzle.geometry.bSpline.knots = knots;
nozzle.geometry.bSpline.coefs = coefs;

% ============= PRINT DESIGN VARIABLES FOR TROUBLESHOOTING ===============

fprintf('%0.4f ',coefs(1,:));
fprintf('\n');
fprintf('%0.4f ',coefs(2,:));
fprintf('\n');

% ====================== SET CONSTRAINT CONSTANTS ========================
minRequiredThrust = 25000; % N

% ========================== PREPARE INPUTS ==============================

% ------------------------ SET FLIGHT REGIME -----------------------------
% High speed, high altitude case
altitude = 35000; % in feet
mach = 0.9;
nozzle.inlet.Tstag = 1021.5;
nozzle.inlet.Pstag = 1.44925e5;

% --------------------- SET HEAT TRANSFER PARAMS -------------------------
nozzle.hInf = 25; % W/m^2-K, heat transfer coefficient from external nozzle wall to environment

% --------------------- SET MATERIAL PROPERTIES --------------------------
% SEPCARBINOX A500, a ceramic matrix composite
nozzle.wall.k = 8.6; % W/m*K, thermal conductivity of nozzle wall
nozzle.wall.coeffThermalExpansion = 2.3e-6; % 1/K, coefficient of thermal expansion of nozzle wall
nozzle.wall.E = 80e9; % Pa, elastic modulus of nozzle wall
nozzle.wall.poissonRatio = 0.3; % Poisson ratio of nozzle wall

% ------------------- SET err TOLERANCE RANGES -------------------------
err.betweenIterations.inletMach = 1e-10;
err.solver.inletMach = 1e-8;
err.betweenIterations.exitTemp = 1e-6;
err.solver.apparentThroatLocation = 1e-6;
err.solver.M2relative = 1e-10;
err.solver.M2absolute = 1e-10;
err.dMdxDenominator = 4; % this is not an err tolerance, rather it is used to set the slope of dMdx in the transonic regime

% ----------------------- SET NOZZLE GEOMETRY ----------------------------

nozzle.geometry.shape = 'B-spline-mex';

% knots and coefs defined at beginning of file

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

% -------------------- SET NOZZLE WALL GEOMETRY --------------------------
nozzle.wall.shape = 'piecewise-linear';
% nozzle.wall.seed = an array of the form [x; y] where [x,y]
% denote the location of the control points with the origin being at
% the center of the inlet area
% nozzle.wall.breaks = vector giving location of breaks in piecewise
% function
nozzle.wall.seed = [0, 0.01; 
                    nozzle.geometry.xThroat, 0.01; 
                    nozzle.geometry.length, 0.01];
nozzle.wall.breaks = [0;
                      nozzle.geometry.xThroat;
                      nozzle.geometry.length];
                  
% --------------------- FLUID INPUT PROPERTIES ---------------------------
fluid.gam = 1.4; % ratio of specific heats
fluid.R = 287.06; % J/kg-K, specific gas constant

% ==================== RETURN VOLUME IF REQUESTED ========================
if( strcmp(output,'volume') )
    volume = wallVolume(knots,coefs,nozzle.wall,'B-spline');
    varargout = {volume};
    return;
end

% ================= OTHERWISE RUN NOZZLE CALCULATIONS ====================

% ----------------- CALCULATE FREESTREAM PROPERTIES ----------------------
atm = StndAtm(altitude*0.3048,'SI'); % obtain standard atmosphere characteristics
freestream.P = atm.P; % Pa, atmospheric pressure
freestream.T = atm.T; % K, atmospheric temperature
freestream.M = mach;
freestream.U = freestream.M*sqrt(fluid.gam*fluid.R*freestream.T);

if( strcmp(fidelity,'low') )
    [ nozzle ] = nozzleNonIdeal( fluid, freestream, nozzle, err );  
elseif( strcmp(fidelity,'med') )
    nozzle.meshSize  = 'coarse'; % 'coarse', 'medium', or 'fine'
    nozzle.governing = 'euler'; % 'euler' or 'rans'
    
    % Run calculations in a separate directory
    
    % 1) Add nozzle code to Matlab's path
    addpath(F100dir);
    % 2) Determine where this optimization code is running
    homedir = pwd;
    % 3) Create a new subdirectory to run SU2 calculations in & move there
    dirname = tempname; % Matlab generates a temp name in the OS's temp directory (if one exists)
    indexUniqueID = regexp(dirname,tempdir,'end') + 1; % find unique ID portion of dirname
    dirname = dirname(indexUniqueID:end); % remove OS temp folder's name from beginning
    mkdir(dirname);
    cd(dirname);
    workdir = pwd; % where SU2 will work
    % 4) Run SU2
    try
        [ nozzle ] = nozzleCFD( fluid, freestream, nozzle, err );
    catch
        error('nozzleCFD did not complete successfully');
    end
    % 5) Remove files generated by SU2
    files = strsplit(ls);
    for ii = 1:length(files)
        if(strcmp(files{ii},'')) % ignore empty strings
            continue;
        else
            delete(files{ii});
        end
    end
    % 6) Remove subdirectory
    cd('..');
    rmdir(dirname);  
else
    error('fidelity not supported');
end

% ----------------- PRINT SOME USEFUL DATA TO SCREEN ---------------------
fprintf('Nozzle Volume: %0.4f cm^3\n',nozzle.geometry.volume*100^3);
fprintf('Nozzle Est. Thrust: %0.4f N\n',nozzle.netThrust);
fprintf('Nozzle Mass Flow Rate: %0.4f kg/s\n',nozzle.massFlowRate);

% ============================ OUTPUT DATA ===============================

if (strcmp(output,'nonlcon')) % nonlinear constraint functions for this case
    Ceq = [];
    
    C(1) = minRequiredThrust - nozzle.netThrust; % should be <= 0
    
    varargout = {C Ceq}; % required return format by fmincon

else % output everything
    varargout = {nozzle};
end


