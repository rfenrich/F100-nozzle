% Test nozzleGeometry function by comparing areas calculated from the A(x)
% area function, D(x) diameter function, and dAdx(x) area change function.
%
% Rick Fenrich 7/3/15

% Nozzle properties
nozzleLength = 1.2;
xThroat = 0.25;
ARatio1 = 5;
ARatio2 = 4;
Dinlet = 0.88;

% Make necessary functions
A = @(x) nozzleGeometry(x,Dinlet,nozzleLength,xThroat,ARatio1,ARatio2,'linear','A');
dAdx = @(x) nozzleGeometry(x,Dinlet,nozzleLength,xThroat,ARatio1,ARatio2,'linear','dAdx');
D = @(x) nozzleGeometry(x,Dinlet,nozzleLength,xThroat,ARatio1,ARatio2,'linear','D');
t = @(x) 0.01; % m, thickness of wall

% Discretize x-axis
xPosition = linspace(0,nozzleLength,100);

% Calculate areas based on diameter D, and dAdx
AfromD = pi*D(xPosition).^2/4;
AfromdAdx = pi*Dinlet^2/4 + cumtrapz(xPosition,dAdx(xPosition));

% Plot results - all calculated areas should be the same
figure; hold on;
plot(xPosition, A(xPosition),'b-');
plot(xPosition, AfromD,'r:');
plot(xPosition, AfromdAdx,'g-.');
legend('A','A from D','A from dAdx');
    