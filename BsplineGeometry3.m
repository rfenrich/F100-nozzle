function [ varargout ] = BsplineGeometry3( x, outputResult, knots, coefs )
% BsplineGeometry3.m calculates the nozzle radius, diameter, area, change
% in area, and thickness for a 3rd degree NURB spline (B-spline) 
% parameterized nozzle geometry. The knots vector and coefficients matrix 
% are used to define the B-spline, and have several requirements (see 
% below under general usage). Use consistent units.
%
% For faster calculations, use BsplineGeometry3Mex instead. To do this, 
% compile BsplineGeometry3Mex.cpp in Matlab with the command:
% >> mex BsplineGeometry3Mex.cpp
%
% General usage is:
% [output] = BsplineGeometry3(x, outputResult, knots, coefs) where:
% x = dimensional location along nozzle length (scalar or vector)
% outputResult = string noting which result is desired to be output (can 
%       take the value of 'A', 'dAdx', 'y', 'D', or 't'. If
%       outputResult is empty, then A, dAdx, D, and t are all output.
% knots = knots vector of B-spline. The elements of the knots vector should
%       be increasing integers, spaced by either 0 or 1 from their 
%       neighbors. The first 3 knots (elements) and the last 4 knots should
%       have the same value.
% coefs = coefficient matrix of B-spline. Should be a [2 x C] matrix where
%       C is the number of coefficients. The first row in the coefficients
%       matrix should contain the x-coordinates of the control points, and
%       the second row should contain the y-coordinates of the control
%       points. Repeat the first and last coefficients once.
%
% To ensure a 3rd degree NURB spline is created, C = K - P - 1, where C is
% the number of control points (coefficients), K is the number of knots,
% and P is the order, in this case P = 3. Since 3 knots are repeated at the
% beginning and 4 knots at the end of the knots vector, the NURB spline 
% will begin at the first control point and end at the last control point.
%
% Example usage:
% >> knots = [0 0 0 1 2 3 4 4 4 4]';
% >> coefs = [0 0 0.05  0.27  0.65  1       1;
%             0 0 0    -0.06 -0.02  0.0038  0.0038];
% >> diameter = BsplineGeometry3([0 0.25 0.75 1], 'D', knots, coefs);
% 
% Furthermore, the function will provide the nozzle throat location(s) and
% size(s) when called in the following way:
% [xThroat, yThroat] = BsplineGeometry3(x, 'throat', knots, coefs)
%
% OUTPUTS:
% A = area corresponding to nozzle cross-sections at location x
% dAdx = rate of change of area corresponding to " "
% y = radius corresponding to " "
% D = diameter corresponding to " "
% xThroat = x-coordinate of throat
% yThroat = y-coordinate of throat
%
% Rick Fenrich 2/8/16

% Properties of B-spline
c = length(coefs); % number of control points
k = length(knots); % knots
p = k-c-1; % spline degree

if(p ~= 3)
    error('Calculations are only for a spline of degree 3');
end

% Determine x value at breaks
uniqueKnots = knots(3:end-3); % for repeated knots at start and end
[xKnot, ~, ~, ~] = uMap(knots,coefs,uniqueKnots);

% At first, lump all results together (i.e. calculate everything)
tolerance = 1e-6; % tolerance for Newton solver

% Determine y value at every x value
y = zeros(length(x),1);
dydx = zeros(length(x),1);
xKnot(end) = xKnot(end) + 1e-6; % so the following algorithm works at last control point
for ii = 1:length(x)

    % Determine lower and upper bounds on u
    s = find(x(ii)<xKnot,1,'first') - 1; % segment number
    uLower = knots(s + p - 1);
    uUpper = knots(s + p);
    
    % Pick a guess for u
    u = (x(ii) - xKnot(s))/(xKnot(s+1) - xKnot(s))*(uUpper - uLower) + uLower;
    %u = (uLower + uUpper)/2;
    
    % Calculate x and dxdu corresponding to u
    [xEst, ~, dxduEst, ~] = uMap(knots,coefs,u);
    
    % Perform 1 Newton iteration
    if(dxduEst == 0);
        uNew = 0;
    else
        uNew = u - (xEst - x(ii))/dxduEst;
    end
    
    % Perform remaining Newton iterations
    counter = 0;
    while( abs((uNew-u)/uNew) > tolerance)
        
        u = uNew;

        [xEst, ~, dxduEst, ~] = uMap(knots,coefs,u);
				
        if(dxduEst == 0)
            uNew = u;
        else
            uNew = u - (xEst - x(ii))/dxduEst;
        end

        counter = counter + 1;
        if(counter > 10)
            break;
        end
    end
    
    u = uNew;
    [~, yEst, dxduEst, dyduEst] = uMap(knots,coefs,u);

    % Calculate y corresponding to given x
    y(ii) = yEst;

    % Calculate dydx
    if(dxduEst == 0)
        dydx(ii) = 0;
    else
        dydx(ii) = dyduEst/dxduEst;
    end

end

% ========================== OUTPUT RESULTS ==============================

if(strcmp(outputResult,'A'))
    
    varargout{1} = pi*y.^2;
    
elseif(strcmp(outputResult,'dAdx'))
    
    varargout{1} = 2*pi*y.*dydx;
    
elseif(strcmp(outputResult,'y'))
    
    varargout{1} = y;
    
elseif(strcmp(outputResult,'D'))
    
    varargout{1} = 2*y;    
    
elseif(strcmp(outputResult,'t'))
    % For right now, set a constant thickness for the nozzle, regardless of
    % shape. Can incorporate variable thickness calculations later above:
    
    varargout{1} = 0.01*ones(length(x),1);
    
elseif(strcmp(outputResult,'throat')) % find location & height of throat
    
    xSegMin = zeros(length(xKnot)-1,1);
    ySegMin = zeros(length(xKnot)-1,1);
    for ii = 1:length(xKnot)-1 % for each segment
        throatFunc = @(r) utoyMap(knots,coefs,r);
        uSegMin = fminbnd(throatFunc,uniqueKnots(ii),uniqueKnots(ii+1));
        xSegMin(ii) = utoxMap(knots,coefs,uSegMin);
        ySegMin(ii) = utoyMap(knots,coefs,uSegMin);
    end
    
    [val,ind] = min(ySegMin);
    xThroat = xSegMin(ind);
    yThroat = val;
    
    varargout{1} = xThroat;
    varargout{2} = yThroat; % throat diameter
    
else
    
    varargout{1} = pi*y.^2; % Area
    varargout{2} = 2*pi*y.*dydx; % dAdx
    varargout{3} = 2*y; % D
    varargout{4} = 0.01*ones(length(x),1); % thickness, constant for now
    
end

end

function [xVector, yVector, dxduVector, dyduVector] = uMap(knots,coefs,uVector)

    xVector = zeros(length(uVector),1);
    yVector = zeros(length(uVector),1);
    dxduVector = zeros(length(uVector),1);
    dyduVector = zeros(length(uVector),1);

    % Fix the end indices of the knot vector, so the highest knot index search
    % works out
    endIndex = find(knots==knots(end));
    knots(endIndex) = knots(endIndex) + 1e-6;

    for qq = 1:length(uVector)

        u = uVector(qq); % given u
 
        % calculate the highest knot index that gives a value of u below the given value (1-based)
        i = find(knots<=u,1,'last');
        
        % nn = min(i,4);
        if(i == 1), nn = 1;
        elseif(i == 2), nn = 2;
        elseif(i == 3), nn = 3;
        else nn = 4;
        end

        point = zeros(2,1);
        pointDeriv = zeros(2,1);
        ii = 0;
        while ii > -nn % for each contributing basis to the calculation

            j = i + ii;

            % redefine k1 through k5 here
            k1 = knots(j);
            k2 = knots(j+1);
            k3 = knots(j+2);
            k4 = knots(j+3);
            k5 = knots(j+4);
            if(ii == 0)
                % calculate basis N1
                if(k1 == k2), Ncurrent = 0;
                else Ncurrent = (u-k1)/(k4-k1)*(u-k1)/(k3-k1)*(u-k1)/(k2-k1);
                    dNducurrent = -(3*(k1 - u)^2)/((k1 - k2)*(k1 - k3)*(k1 - k4));
                end
            elseif(ii == -1)
                % calculate basis N2
                if(k2 == k3), Ncurrent = 0;
                else Ncurrent = (u-k1)/(k4-k1)*((u-k1)/(k3-k1)*(k3-u)/(k3-k2) + (k4-u)/(k4-k2)*(u-k2)/(k3-k2)) + (k5-u)/(k5-k2)*(u-k2)/(k4-k2)*(u-k2)/(k3-k2);
                    dNducurrent = (((k1 - u)*(k3 - u))/((k1 - k3)*(k2 - k3)) + ((k2 - u)*(k4 - u))/((k2 - k3)*(k2 - k4)))/(k1 - k4) + (k2 - u)^2/((k2 - k3)*(k2 - k4)*(k2 - k5)) + ((k5 - u)*(2*k2 - 2*u))/((k2 - k3)*(k2 - k4)*(k2 - k5)) + (2*(k1 - u)*(k1*k2 - k3*k4 - k1*u - k2*u + k3*u + k4*u))/((k1 - k3)*(k1 - k4)*(k2 - k3)*(k2 - k4));
                end
            elseif(ii == -2)
                % calculate basis N3
                if(k3 == k4), Ncurrent = 0;
                else Ncurrent = (u-k1)/(k4-k1)*(k4-u)/(k4-k2)*(k4-u)/(k4-k3) + (k5-u)/(k5-k2)*((u-k2)/(k4-k2)*(k4-u)/(k4-k3) + (k5-u)/(k5-k3)*(u-k3)/(k4-k3));
                    dNducurrent = - (((k2 - u)*(k4 - u))/((k2 - k4)*(k3 - k4)) + ((k3 - u)*(k5 - u))/((k3 - k4)*(k3 - k5)))/(k2 - k5) - (k4 - u)^2/((k1 - k4)*(k2 - k4)*(k3 - k4)) - ((k1 - u)*(2*k4 - 2*u))/((k1 - k4)*(k2 - k4)*(k3 - k4)) - (2*(k5 - u)*(k2*k3 - k4*k5 - k2*u - k3*u + k4*u + k5*u))/((k2 - k4)*(k2 - k5)*(k3 - k4)*(k3 - k5));
                end
            else
                % calculate basis N4
                if(k4 == k5), Ncurrent = 0;
                else Ncurrent = (k5-u)/(k5-k2)*(k5-u)/(k5-k3)*(k5-u)/(k5-k4);  
                    dNducurrent = (3*(k5 - u)^2)/((k2 - k5)*(k3 - k5)*(k4 - k5));
                end
            end

            point = point + coefs(:,j)*Ncurrent;
            pointDeriv = pointDeriv + coefs(:,j)*dNducurrent;

            ii = ii - 1;

            if(ii < -3 - 1), break; end

        end

        xVector(qq) = point(1);
        yVector(qq) = point(2);
        dxduVector(qq) = pointDeriv(1);
        dyduVector(qq) = pointDeriv(2);

    end

end

function y = utoyMap(knots,coefs,u)

    [~, y, ~, ~] = uMap(knots,coefs,u);
    
end

function x = utoxMap(knots,coefs,u)

    [x, ~, ~, ~] = uMap(knots,coefs,u);
    
end

