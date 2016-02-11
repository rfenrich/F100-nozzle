% exampleOptDriver.m runs a deterministic optimization of nozzle shape 
% for maximum thrust using Matlab's fmincon. The nozzle shape is
% is parameterized using a 3rd degree B-spline. There are 14 design
% variables, which correspond to the x and y coordinates of 7 control
% points. Inlet size of the nozzle is fixed. The optimization 
% formulation is:
%
% min                    weight
% s.t.             thrust >= 25,000 N
%       44 kg/s <= mass flow rate <= 46 kg/s  --> for a well posed problem
%                       Ax <= b               --> constraints on relative control point movement
%                    lb <= x <= b             --> constraints on individual control point movement
%       where x is the vector of design variables
%
% Rick Fenrich 2/11/16

clear all; close all; clc;

% Set up initial B-spline shape
Dinlet = 0.3255*2;
knots = [0 0 0 0 1 2 3 4 5 6 7 8 9 9 9 9]';
coefs = [0.0000 0.0000 0.1470 0.1578 0.1707 0.2180 0.2233 0.3186 0.4178 0.6008 0.6700 0.6700;
         0.3255 0.3255 0.3255 0.3255 0.3255 0.3250 0.2945 0.2719 0.2944 0.3056 0.3050 0.3050];

% Set up design variables (x and y coordinates of middle 2 control points, 
% and y coordinate of last control point)
x(1:7) = coefs(1,6:12);
x(8:14) = coefs(2,6:12);

% Set up ranges of design variables
lb = 0.8*x;
ub = 1.2*x;

n = length(x); % number of design variables

% Set up linear inequality constraints
A = zeros(11,length(x)); % 11 constraints total
b = zeros(11,1); 
delta = 1e-3; % delta for spacing between control points, so code doesn't break

% x and y position constraints
A(1,1) = -1; b(1) = -coefs(1,5) - delta*10; % 6th control point should remain to right of 5th
A(2,1) = 1; A(2,2) = -1; b(2) = -delta*10; % 7th control point should remain to right of 6th c.p.
A(3,2) = 1; A(3,3) = -1; b(3) = -delta*10; % 8th c.p.  " "
A(4,3) = 1; A(4,4) = -1; b(4) = -delta*10; % 9th c.p.  " "
A(5,4) = 1; A(5,5) = -1; b(5) = -delta*10; % 10th c.p.  " "
A(6,5) = 1; A(6,6) = -1; b(6) = -delta*10; % 11th c.p.  " "

A(7,10) = 1; A(7,9) = -1; b(7) = -delta; % 8th c.p. should remain below 7th c.p.
A(8,9) = 1; b(8) = Dinlet/2; % 7th c.p. should be below inlet height
A(9,10) = 1; A(9,8) = -1; b(9) = -delta; % 8th c.p. should remain below 6th c.p.
A(10,8) = 1; b(10) = Dinlet/2; % 6th c.p. should be below inlet height

A(11,10) = 1; A(11,11) = -1; b(11) = -delta; % 9th c.p. should remain above 8th c.p.
A(12,10) = 1; A(12,12) = -1; b(12) = -delta*10; % 10th c.p. should remain above 8th c.p.
A(13,10) = 1; A(13,13) = -1; b(13) = -delta*10; % 11th c.p. should remain above 8th c.p.

% Set linear equality constraints so that nozzle inlet is fixed size
Aeq = zeros(1,length(x));
Aeq(1,6) = 1; Aeq(1,7) = -1; beq(1) = 0; % 11th and 12th control point x-coordinate is the same
Aeq(2,13) = 1; Aeq(2,14) = -1; beq(2) = 0; % 11th and 12th control point y-coordinate is the same

% Set nonlinear inequality constraint function
nonlconFun = @(r) exampleWrapper2(r,knots,coefs,'nonlcon');

% Set optimization options
options.MaxIter = 50;
options.Algorithm = 'sqp';
options.Display = 'iter-detailed';
options.FinDiffRelStep = 1e-3;
opttions.TolX = 1e-10;

% Print data to screen
fprintf('Number design variables: %i\n',n);

objFun = @(r) exampleWrapper2(r,knots,coefs,'volume');

figure;
hold on;
tic;
[sol,val] = fmincon(objFun,x,A,b,Aeq,beq,lb,ub,nonlconFun,options);
timeToEnd = toc;
fprintf('Time to completion: %0.2f sec\n',timeToEnd);


