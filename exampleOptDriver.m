% exampleOptDriver.m runs a deterministic optimization of nozzle shape 
% for maximum thrust using Matlab's fmincon. The nozzle shape is
% is parameterized using a NURB spline (B-spline). There are 9 design
% variables, which correspond to the x and y coordinates of 5 control
% points. Length and inlet size of the nozzle is fixed. The optimization 
% formulation is:
%
% min                    weight
% s.t.             thrust >= 25,000 N
%       44 kg/s <= mass flow rate <= 46 kg/s  --> for a well posed problem
%                       Ax <= b               --> constraints on relative control point movement
%                    lb <= x <= b             --> constraints on individual control point movement
%       where x is the vector of design variables
%
% Rick Fenrich 2/4/16

clear all; close all; clc;

% Set up initial B-spline shape
Dinlet = 0.3255*2;
knots = [0 0 0 1 2 3 4 5 5 5]';
coefs = [0         0.1720    0.2095    0.2328    0.3170    0.4997   0.6700;
         0.3255    0.3255    0.3254    0.2935    0.2733    0.3044   0.3050];

% Set up design variables (x and y coordinates of middle 2 control points, 
% and y coordinate of last control point)
x(1:4) = coefs(1,3:6);
x(5:9) = coefs(2,3:7);

% Set up ranges of design variables
lb = 0.8*x;
ub = 1.2*x;

n = length(x); % number of design variables

% Set up linear inequality constraints
A = zeros(11,9); % 11 constraints total
b = zeros(11,1); 
delta = 1e-3; % delta for spacing between control points, so code doesn't break

% x and y position constraints
A(1,1) = -1; b(1) = -coefs(1,2)*10; % 3rd control point should remain to right of 2nd
A(2,1) = 1; A(2,2) = -1; b(2) = -delta*10; % 4th control point should remain to right of 3rd c.p.
A(3,2) = 1; A(3,3) = -1; b(3) = -delta*10; % 5th c.p.  " "
A(4,3) = 1; A(4,4) = -1; b(4) = -delta*10; % 6th c.p.  " "
A(5,4) = 1; b(5) = coefs(1,end)-delta*10; % 6th c.p. should remain to left of 7th c.p.
A(6,7) = 1; A(6,6) = -1; b(6) = -delta; % 5th c.p. should remain below 4th c.p.
A(7,6) = 1; b(7) = Dinlet/2; % 4th c.p. should be below inlet height
A(8,7) = 1; A(8,5) = -1; b(8) = -delta; % 5th c.p. should remain below 3rd c.p.
A(9,5) = 1; b(9) = Dinlet/2; % 3rd c.p. should be below inlet height
A(10,7) = 1; A(10,8) = -1; b(10) = -delta; % 6th c.p. should remain above 5th c.p.
A(11,7) = 1; A(11,9) = -1; b(11) = -delta*10; % 7th c.p. should remain above 5th c.p.

% Set linear equality constraints so that nozzle inlet is fixed size
beq = [];
Aeq = [];

% Set nonlinear inequality constraint function
nonlconFun = @(r) exampleWrapper(r,knots,coefs,'nonlcon');

% Set optimization options
options.MaxIter = 50;
options.Algorithm = 'sqp';
options.Display = 'iter-detailed';
options.FinDiffRelStep = 1e-3;
opttions.TolX = 1e-10;
%options.FunValCheck = 'on';

% Print data to screen
fprintf('Number design variables: %i\n',n);

objFun = @(r) exampleWrapper(r,knots,coefs,'volume');

figure;
hold on;
tic;
[sol,val] = fmincon(objFun,x,A,b,Aeq,beq,lb,ub,nonlconFun,options);
timeToEnd = toc;
fprintf('Time to completion: %0.2f sec\n',timeToEnd);


