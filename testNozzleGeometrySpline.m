% Test splineGeometry function by comparing areas calculated from the A(x)
% area function, D(x) diameter function, and dAdx(x) area change function.
%
% Rick Fenrich 9/1/15

% Nozzle properties
nozzleLength = 1.2;
xThroat = 0.25;
ARatio1 = 5;
ARatio2 = 4;
Dinlet = 0.88;

% Create nozzle geometry from which splines should be based:
A = @(x) nozzleGeometry(x,'A',Dinlet,nozzleLength,xThroat,ARatio1,ARatio2,'linear');
dAdx = @(x) nozzleGeometry(x,'dAdx',Dinlet,nozzleLength,xThroat,ARatio1,ARatio2,'linear');
D = @(x) nozzleGeometry(x,'D',Dinlet,nozzleLength,xThroat,ARatio1,ARatio2,'linear');
t = @(x) 0.01; % m, thickness of wall

% Set control points for splines (xNode), value of function at control
% point (yNode), and slopes at start and end of spline (slopes)
xNode = linspace(0,nozzleLength,50);
yNode = D(xNode)/2;
slopes = [0, 0];
pp = spline(xNode,[slopes(1); yNode; slopes(2)]); % perform piecewise cubic spline interpolation

% Make necessary functions for splined nozzle shape
A = @(x) splineGeometry(x,'A',pp);
dAdx = @(x) splineGeometry(x,'dAdx',pp);
D = @(x) splineGeometry(x,'D',pp);
t = @(x) 0.01; % m, thickness of wall

% Discretize x-axis
xPosition = linspace(0,nozzleLength,200);

% Calculate areas based on diameter D, and dAdx
AfromD = pi*D(xPosition).^2/4;
AfromdAdx = pi*Dinlet^2/4 + cumtrapz(xPosition,dAdx(xPosition));

% Plot results - all calculated areas should be the same
figure; hold on;
plot(xPosition, A(xPosition),'b-');
plot(xPosition, AfromD,'r:');
plot(xPosition, AfromdAdx,'g-.');
legend('A','A from D','A from dAdx');
    