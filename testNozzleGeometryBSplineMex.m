% Test BsplineGeometry function by comparing areas calculated from the A(x)
% area function, D(x) diameter function, and dAdx(x) area change function.
%
% Rick Fenrich 1/19/16

clear all; close all; clc;

% Nozzle properties
nozzleLength = 1;
xThroat = 0.3;
ARatio1 = 1.368;
ARatio2 = 1.4;
Dinlet = 0.651;

% Create linear nozzle geometry from which splines should be based:
DLinear = @(x) nozzleGeometry(x,'D',Dinlet,nozzleLength,xThroat,ARatio1,ARatio2,'linear');

% Create NURBS spline geometry
knots = [0 0 0 1 2 3 4 4 4]';
coefs = [0 0 0.05 0.27  0.65  1;
         0 0 0   -0.06   -0.02     0.0038];
coefs(1,:) = coefs(1,:)*nozzleLength;
coefs(2,:) = coefs(2,:)+Dinlet/2;

% Set bounds on coefficients
lb = coefs(2,:) - 0.8*(coefs(2,:) + Dinlet/2);
ub = coefs(2,:) + 0.8*(coefs(2,:) + Dinlet/2);

% Make necessary functions for splined nozzle shape
A = @(x) BsplineGeometryMex(x,1,knots,coefs');
dAdx = @(x) BsplineGeometryMex(x,2,knots,coefs');
D = @(x) BsplineGeometryMex(x,3,knots,coefs');
t = @(x) 0.01; % m, thickness of wall

% Discretize x-axis
xPosition = linspace(0,nozzleLength,200)';

% Calculate areas based on diameter D, and dAdx
AfromD = pi*D(xPosition).^2/4;
AfromdAdx = pi*Dinlet^2/4 + cumtrapz(xPosition,dAdx(xPosition));

% Plot shape of nozzle
hold on;
plot(xPosition,DLinear(xPosition)/2);
plot(xPosition,D(xPosition)/2);
plot(coefs(1,:),coefs(2,:),'ko');
axis equal

% Plot results - all calculated areas should be the same
figure; hold on;
plot(xPosition, A(xPosition),'b-');
plot(xPosition, AfromD,'r:');
plot(xPosition, AfromdAdx,'g-.');
legend('A','A from D','A from dAdx');

    