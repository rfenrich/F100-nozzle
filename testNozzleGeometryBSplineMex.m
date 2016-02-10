% Test BsplineGeometry function by comparing areas calculated from the A(x)
% area function, D(x) diameter function, and dAdx(x) area change function.
%
% Rick Fenrich 2/10/16

clear all; close all; clc;

% Nozzle properties
nozzleLength = 1;
Dinlet = 0.651;

%% 2nd degree B-spline

% Create NURBS spline geometry
knots = [0 0 0 1 2 3 4 4 4]';
coefs = [0 0 0.05 0.27  0.65  1;
         0 0 0   -0.06   -0.02     0.0038];
coefs(1,:) = coefs(1,:)*nozzleLength;
coefs(2,:) = coefs(2,:)+Dinlet/2;

% Make necessary functions for splined nozzle shape
A = @(x) BsplineGeometryMex(x,1,knots,coefs');
dAdx = @(x) BsplineGeometryMex(x,2,knots,coefs');
D = @(x) BsplineGeometryMex(x,3,knots,coefs');
t = @(x) 0.01; % m, thickness of wall

x = linspace(0,1,100)';
Dout = D(x);
Aout = A(x);
dAdxout = dAdx(x);

% Discretize x-axis
xPosition = linspace(0,nozzleLength,200)';

% Calculate areas based on diameter D, and dAdx
AfromD = pi*D(xPosition).^2/4;
AfromdAdx = pi*Dinlet^2/4 + cumtrapz(xPosition,dAdx(xPosition));

% Plot shape of nozzle
figure; hold on;
f = spmak(knots,coefs); % compare with Matlab
fnplt(f);
plot(xPosition,D(xPosition)/2);
plot(coefs(1,:),coefs(2,:),'ko');
axis equal
legend('Matlab','self');
title('2nd degree spline');

% Plot results - all calculated areas should be the same
figure; hold on;
plot(xPosition, A(xPosition),'b-');
plot(xPosition, AfromD,'r:');
plot(xPosition, AfromdAdx,'g-.');
legend('A','A from D','A from dAdx');
title('2nd degree spline');

%% 3rd degree B-spline

% Nozzle properties
nozzleLength = 0.67;
Dinlet = 0.651;

% Create NURBS spline geometry
knots = [0 0 0 0 1 2 3 4 4 4 4]';
coefs = [0         0         0.2095    0.2328    0.3170    0.6700   0.6700;
         0.3255    0.3255    0.3254    0.2935    0.2733    0.3050   0.3050];

% Make necessary functions for splined nozzle shape
A = @(x) BsplineGeometryMex(x,1,knots,coefs');
dAdx = @(x) BsplineGeometryMex(x,2,knots,coefs');
D = @(x) BsplineGeometryMex(x,3,knots,coefs');
t = @(x) 0.01; % m, thickness of wall

x = linspace(0,1,100)';
Dout = D(x);
Aout = A(x);
dAdxout = dAdx(x);

% Discretize x-axis
xPosition = linspace(0,nozzleLength,200)';

% Calculate areas based on diameter D, and dAdx
AfromD = pi*D(xPosition).^2/4;
AfromdAdx = pi*Dinlet^2/4 + cumtrapz(xPosition,dAdx(xPosition));

% Plot shape of nozzle
figure; hold on;
f = spmak(knots,coefs); % compare with Matlab
fnplt(f);
plot(xPosition,D(xPosition)/2);
plot(coefs(1,:),coefs(2,:),'ko');
legend('Matlab','self');
axis equal
title('3rd degree spline');

% Plot results - all calculated areas should be the same
figure; hold on;
plot(xPosition, A(xPosition),'b-');
plot(xPosition, AfromD,'r:');
plot(xPosition, AfromdAdx,'g-.');
legend('A','A from D','A from dAdx');
title('3rd degree spline');



    