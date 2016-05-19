% exampleOptDriver.m runs a deterministic optimization of nozzle shape 
% for minimum weight using Matlab's fmincon. The nozzle shape is
% is parameterized using a 3rd degree B-spline. There are 21 design
% variables, each of which corresponds to an x or y coordinate of a control
% point. The inlet size of the nozzle is fixed. The optimization 
% formulation is:
%
% min                    weight
% s.t.             thrust >= 25,000 N
%                       Ax <= b
%                    lb <= x <= b
%       where x is the vector of design variables
%
% Note that no mass flow rate constraint is necessary. The linear
% constraints Ax <= b include constraints which:
% a) prohibit control points from crossing over each other in the
% x-direction
% b) require a converging pre-throat area
% c) require the throat to be the lowest point
% d) slope of line drawn between two adjacent control points to not be too
% steep
%
% Results using low-fidelity nozzle:
% Volume at optimum: 0.009844 m^3
% Thrust at optimum: 25,000 N
% Time to completion: 796 sec on 4 cores
% Number of iterations: 20
%
% Rick Fenrich 3/17/16 updated 5/17/16

clear all; close all; clc;

fidelity = 'med'; % 'low' or 'med'
F100dir = '/home/rick/Documents/Research/F100-nozzle'; % directory of F100-nozzle code

% Set up initial B-spline shape
Dinlet = 0.3255*2;
knots = [0 0 0 0 1:14 15 15 15 15]';
coefs = [0.0000 0.0000 0.1500 0.1700 0.1900 0.2124 0.2269 0.2734 0.3218 0.3218 0.3230 0.3343 0.3474 0.4392 0.4828 0.5673 0.6700 0.6700;
         0.3255 0.3255 0.3255 0.3255 0.3255 0.3238 0.2981 0.2817 0.2787 0.2787 0.2787 0.2797 0.2807 0.2936 0.2978 0.3049 0.3048 0.3048];

% Set up design variables
x = zeros(21,1);
x(1:11) = [coefs(1,6:9), coefs(1,11:17)]; % x-coordinates
x(12:21) = [coefs(2,6:9), coefs(2,12:17)]; % y-coordinates

% Set up ranges of design variables
lb = 0.8*x;
ub = 1.2*x;

n = length(x); % number of design variables

% Define linear inequality constraints, Ax <= b
A = zeros(44,length(x)); % 44 constraints total
b = zeros(44,1);
Rinlet = 0.3255;
Xinlet = 0.19;
delta = 1e-3;
delta2 = 1e-3; % for proximity to throat
mMin = -1.8;
mMax = 1.8;
nLinearCon = 44;

conNum = 1;

% Set constraints on x-position
A(1,1) = -1; b(1) = -Xinlet - delta; conNum = conNum + 1;
for ii = conNum:conNum+9
  A(ii,ii) = -1; A(ii,ii-1) = 1; b(ii) = -delta; conNum = conNum + 1;
end

% Set constraints on convergence of pre-throat area
jj = 12;
for ii = conNum:conNum+3
  A(ii,jj) = 1; b(ii) = Rinlet; conNum = conNum + 1; jj = jj + 1;
end

% Set constraints making throat the lowest point
jj = 12;
for ii = conNum:conNum+2
  A(ii,15) = 1; A(ii,jj) = -1; b(ii) = -delta2; conNum = conNum + 1; jj = jj + 1;
end
jj = 16;
for ii = conNum:conNum+5
  A(ii,15) = 1; A(ii,jj) = -1; b(ii) = -delta2; conNum = conNum + 1; jj = jj + 1;
end

% Set constraints for steepness of slopes
A(conNum,1) = mMin; A(conNum,12) = -1; b(conNum) = mMin*Xinlet - Rinlet; conNum = conNum + 1;
A(conNum,12) = 1; A(conNum,1) = -mMax; b(conNum) = Rinlet - mMax*Xinlet; conNum = conNum + 1;
jj = 2;
for ii = conNum:conNum+2
  A(ii,jj) = mMin; A(ii,jj-1) = -mMin; A(ii,jj+11) = -1; A(ii,jj+10) = 1; conNum = conNum + 1; jj = jj + 1;
end
jj = 2;
for ii = conNum:conNum+2
  A(ii,jj+11) = 1; A(ii,jj+10) = -1; A(ii,jj) = -mMax; A(ii,jj-1) = mMax; conNum = conNum + 1; jj = jj + 1;
end
jj = 6;
for ii = conNum:conNum+5
  A(ii,jj) = mMin; A(ii,jj-1) = -mMin; A(ii,jj+10) = -1; A(ii,jj+9) = 1; conNum = conNum + 1; jj = jj + 1;
end
jj = 6;
for ii = conNum:conNum+5
  A(ii,jj+10) = 1; A(ii,jj+9) = -1; A(ii,jj) = -mMax; A(ii,jj-1) = mMax; conNum = conNum + 1; jj = jj + 1;
end

% Check feasibility of linear constraints at starting point
b2 = A*x;
bDiff = b - b2;
nConInfeasible = 0;
for ii = 1:length(bDiff)
    if( bDiff(ii) < 0 )
        nConInfeasible = nConInfeasible + 1;
    end
end
if( nConInfeasible > 0 )
    fprintf('%i linear inequality constraints infeasible at starting point\n',nConInfeasible);
else
    fprintf('no linear inequality infeasibilities!\n');
end

% Set nonlinear inequality constraint function
nonlconFun = @(r) exampleWrapper(r,knots,coefs,fidelity,'nonlcon',F100dir);

% Set optimization options
options.MaxIter = 50;
options.Algorithm = 'sqp';
options.Display = 'iter-detailed';
options.FinDiffRelStep = 1e-3;
opttions.TolX = 1e-6;
options.UseParallel = true; % estimate finite difference gradients in parallel

% Print data to screen
fprintf('Number design variables: %i\n',n);

objFun = @(r) exampleWrapper(r,knots,coefs,fidelity,'volume',F100dir);

figure;
hold on;
tic;
[sol,val] = fmincon(objFun,x,A,b,[],[],lb,ub,nonlconFun,options);
timeToEnd = toc;
fprintf('Time to completion: %0.2f sec\n',timeToEnd);

% nozzle = exampleWrapper(x,knots,coefs,fidelity,'all',F100dir);


