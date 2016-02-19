% exampleOptDriver.m runs a deterministic optimization of nozzle shape 
% for maximum thrust using Matlab's fmincon. The nozzle shape is
% is parameterized using a 3rd degree B-spline. There are 22 design
% variables, which correspond to the x and y coordinates of 11 control
% points. Inlet size of the nozzle is fixed. The optimization 
% formulation is:
%
% min                    weight
% s.t.             thrust >= 25,000 N
%       30 kg/s <= mass flow rate <= 46 kg/s  --> for a well posed problem
%                  min(wall slope) >= -1.2
%                  max(wall slope) <= 1.2
%                       Ax <= b               --> constraints on relative control point movement
%                    lb <= x <= b             --> constraints on individual control point movement
%       where x is the vector of design variables
%
% Converged in 390 seconds & 6 iterations
%
% Rick Fenrich 2/17/16

clear all; close all; clc;

% Set up initial B-spline shape
Dinlet = 0.3255*2;
knots = [0 0 0 0 1:12 13 13 13 13]';
coefs = [0.0000 0.0000 0.1500 0.1700 0.1900 0.2124 0.2269 0.2734 0.3218 0.3343 0.3474 0.4392 0.4828 0.5673 0.6700 0.6700;
         0.3255 0.3255 0.3255 0.3255 0.3255 0.3238 0.2981 0.2817 0.2787 0.2790 0.2804 0.2936 0.2978 0.3049 0.3048 0.3048];

% Set up design variables (x and y coordinates of middle 2 control points, 
% and y coordinate of last control point)
x(1:11) = coefs(1,6:16);
x(12:22) = coefs(2,6:16);

% Set up ranges of design variables
lb = 0.8*x;
ub = 1.2*x;

n = length(x); % number of design variables

% Set up linear inequality constraints
A = zeros(23,length(x)); % 11 constraints total
b = zeros(23,1); 
delta = 1e-2; % delta for spacing between control points, so code doesn't break

% Control points cannot crossover
A(1,1) = -1; b(1) = -coefs(1,5) - delta; % 6th control point should remain to right of 5th
for ii = 2:10
    A(ii,ii) = -1; A(ii,ii-1) = 1; b(ii) = -delta;
end

% Pre-throat area must not diverge
for ii = 11:14
    A(ii,ii+1) = 1; b(ii) = Dinlet/2;
end

% Throat must be lowest point
for ii = 15:17
    A(ii,15) = 1; A(ii,ii-3) = -1; b(ii) = delta;
end
for ii = 18:23
    A(ii,15) = 1; A(ii,ii-2) = -1; b(ii) = delta;
end

% Set linear equality constraints so that nozzle inlet is fixed size
Aeq = zeros(2,length(x));
Aeq(1,10) = 1; Aeq(1,11) = -1; beq(1) = 0; % 15th and 16th control point x-coordinate is the same
Aeq(2,21) = 1; Aeq(2,22) = -1; beq(2) = 0; % 15th and 16th control point y-coordinate is the same

% Set nonlinear inequality constraint function
nonlconFun = @(r) exampleWrapper(r,knots,coefs,'nonlcon');

% Set optimization options
options.MaxIter = 50;
options.Algorithm = 'sqp';
options.Display = 'iter-detailed';
options.FinDiffRelStep = 1e-3;
opttions.TolX = 1e-10;

% Print data to screen
fprintf('Number design variables: %i\n',n);

objFun = @(r) exampleWrapper(r,knots,coefs,'volume');

figure;
hold on;
tic;
[sol,val] = fmincon(objFun,x,A,b,Aeq,beq,lb,ub,nonlconFun,options);
timeToEnd = toc;
fprintf('Time to completion: %0.2f sec\n',timeToEnd);


