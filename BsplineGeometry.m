function [ varargout ] = BsplineGeometry( x, outputResult, knots, coefs )
% BsplineGeometry.m calculates the nozzle radius, diameter, area, change
% in area, and thickness for a 2nd degree NURB spline (B-spline) 
% parameterized nozzle geometry. The knots vector and coefficients matrix 
% are used to define the B-spline, and have several requirements (see 
% below under general usage). Use consistent units.
%
% For faster calculations, use BsplineGeometryMex instead. To do this, 
% compile BsplineGeometryMex.cpp in Matlab with the command:
% >> mex BsplineGeometryMex.cpp
%
% General usage is:
% [output] = BsplineGeometry(x, outputResult, knots, coefs) where:
% x = dimensional location along nozzle length (scalar or vector)
% outputResult = string noting which result is desired to be output (can 
%       take the value of 'A', 'dAdx', 'y', 'D', or 't'. If
%       outputResult is empty, then A, dAdx, D, and t are all output.
% knots = knots vector of B-spline. The elements of the knots vector should
%       be increasing integers, spaced by either 0 or 1 from their 
%       neighbors. The first 3 knots (elements) and the last 3 knots should
%       have the same value.
% coefs = coefficient matrix of B-spline. Should be a [2 x C] matrix where
%       C is the number of coefficients. The first row in the coefficients
%       matrix should contain the x-coordinates of the control points, and
%       the second row should contain the y-coordinates of the control
%       points.
%
% To ensure a 2nd degree NURB spline is created, C = K - P - 1, where C is
% the number of control points (coefficients), K is the number of knots,
% and P is the order, in this case P = 2. Since 3 knots are repeated at the
% beginning and end of the knots vector, the NURB spline will begin at the
% first control point and end at the last control point.
%
% Example usage:
% >> knots = [0 0 0 1 2 3 4 4 4]';
% >> coefs = [0 0 0.05  0.27  0.65  1;
%             0 0 0    -0.06 -0.02  0.0038];
% >> diameter = BsplineGeometry([0 0.25 0.75 1], 'D', knots, coefs);
% 
% Furthermore, the function will provide the nozzle throat location(s) and
% size(s) when called in the following way:
% [xThroat, yThroat] = BsplineGeometry(x, 'throat', knots, coefs)
%
% OUTPUTS:
% A = area corresponding to nozzle cross-sections at location x
% dAdx = rate of change of area corresponding to " "
% y = radius corresponding to " "
% D = diameter corresponding to " "
% xThroat = x-coordinate of throat
% yThroat = y-coordinate of throat
%
% Rick Fenrich 1/29/16

% Properties of B-spline
c = length(coefs); % number of control points
k = length(knots); % knots
p = k-c-1; % spline degree

if(p < 2 || p > 4)
    error('Calculations are only for a spline of degree 2 or 3');
end

% Determine x value at breaks
uniqueKnots = knots(3:end-2); % for 3 repeated knots at start and end
ns = length(uniqueKnots)-1; % number of segments
seg = [1:ns, ns];
xKnot = zeros(length(uniqueKnots),1);
for ii = 1:length(uniqueKnots)
   xKnot(ii) = utoxMap(knots,coefs,seg(ii),uniqueKnots(ii));
end

% ========================== OUTPUT RESULTS ==============================

if(strcmp(outputResult,'A'))
    
    % Determine y value at every x value
    y = zeros(length(x),1);
    xKnot(end) = xKnot(end) + 1e-6; % so the following algorithm works at last control point
    for ii = 1:length(x)

        % Determine segment number for given x
        s = find(x(ii)<xKnot,1,'first') - 1; % segment number   

        % Solve for u using Newton method
        uGuess = s - 0.5; % assuming knots spaced by 1
        for jj = 1:15
           u = uGuess - g(knots,coefs,s,uGuess,x(ii))/dgdu(knots,coefs,s,uGuess);
           if(abs((u-uGuess)/uGuess) <= 1e-6), break; end
           uGuess = u;
        end

        % Calculate y corresponding to given x
        y(ii) = utoyMap(knots,coefs,s,u);

    end

    varargout{1} = pi*y.^2;
    
elseif(strcmp(outputResult,'dAdx'))
    
    % Determine y value at every x value
    y = zeros(length(x),1);
    dydx = zeros(length(x),1);
    xKnot(end) = xKnot(end) + 1e-6; % so the following algorithm works at last control point
    for ii = 1:length(x)

        % Determine segment number for given x
        s = find(x(ii)<xKnot,1,'first') - 1; % segment number   

        % Solve for u using Newton method
        uGuess = s - 0.5; % assuming knots spaced by 1
        for jj = 1:15
           u = uGuess - g(knots,coefs,s,uGuess,x(ii))/dgdu(knots,coefs,s,uGuess);
           if(abs((u-uGuess)/uGuess) <= 1e-6), break; end
           uGuess = u;
        end

        % Calculate y corresponding to given x
        y(ii) = utoyMap(knots,coefs,s,u);

        % Calculate dydx
        dydx(ii) = dydu(knots,coefs,s,u)*dudx(knots,coefs,s,u);

    end

    varargout{1} = 2*pi*y.*dydx;
    
elseif(strcmp(outputResult,'y'))
    
    % Determine y value at every x value
    y = zeros(length(x),1);
    xKnot(end) = xKnot(end) + 1e-6; % so the following algorithm works at last control point
    for ii = 1:length(x)

        % Determine segment number for given x
        s = find(x(ii)<xKnot,1,'first') - 1; % segment number   

        % Solve for u using Newton method
        uGuess = s - 0.5; % assuming knots spaced by 1
        for jj = 1:15
           u = uGuess - g(knots,coefs,s,uGuess,x(ii))/dgdu(knots,coefs,s,uGuess);
           if(abs((u-uGuess)/uGuess) <= 1e-6), break; end
           uGuess = u;
        end

        % Calculate y corresponding to given x
        y(ii) = utoyMap(knots,coefs,s,u);

    end
    
    varargout{1} = y;
    
elseif(strcmp(outputResult,'D'))
    
    % Determine y value at every x value
    y = zeros(length(x),1);
    xKnot(end) = xKnot(end) + 1e-6; % so the following algorithm works at last control point
    for ii = 1:length(x)

        % Determine segment number for given x
        s = find(x(ii)<xKnot,1,'first') - 1; % segment number   

        % Solve for u using Newton method
        uGuess = s - 0.5; % assuming knots spaced by 1
        for jj = 1:15
           u = uGuess - g(knots,coefs,s,uGuess,x(ii))/dgdu(knots,coefs,s,uGuess);
           if(abs((u-uGuess)/uGuess) <= 1e-6), break; end
           uGuess = u;
        end

        % Calculate y corresponding to given x
        y(ii) = utoyMap(knots,coefs,s,u);

    end
    
    varargout{1} = 2*y;    

elseif(strcmp(outputResult,'t'))
    % For right now, set a constant thickness for the nozzle, regardless of
    % shape. Can incorporate variable thickness calculations later above:
    
    varargout{1} = 0.01*ones(length(x),1);
    
elseif(strcmp(outputResult,'throat')) % find location & height of throat
    
    xSegMin = zeros(ns,1);
    ySegMin = zeros(ns,1);
    for ii = 1:ns % for each segment
        throatFunc = @(r) utoyMap(knots,coefs,ii,r);
        uSegMin = fminbnd(throatFunc,uniqueKnots(ii),uniqueKnots(ii+1));
        xSegMin(ii) = utoxMap(knots,coefs,ii,uSegMin);
        ySegMin(ii) = utoyMap(knots,coefs,ii,uSegMin);
    end
    
    [val,ind] = min(ySegMin);
    xThroat = xSegMin(ind);
    yThroat = val;
    
    varargout{1} = xThroat;
    varargout{2} = yThroat; % throat diameter
else
    
    % Determine y value at every x value
    y = zeros(length(x),1);
    dydx = zeros(length(x),1);
    xKnot(end) = xKnot(end) + 1e-6; % so the following algorithm works at last control point
    for ii = 1:length(x)

        % Determine segment number for given x
        s = find(x(ii)<xKnot,1,'first') - 1; % segment number   

        % Solve for u using Newton method
        uGuess = s - 0.5; % assuming knots spaced by 1
        for jj = 1:15
           u = uGuess - g(knots,coefs,s,uGuess,x(ii))/dgdu(knots,coefs,s,uGuess);
           if(abs((u-uGuess)/uGuess) <= 1e-6), break; end
           uGuess = u;
        end

        % Calculate y corresponding to given x
        y(ii) = utoyMap(knots,coefs,s,u);

        % Calculate dydx
        dydx(ii) = dydu(knots,coefs,s,u)*dudx(knots,coefs,s,u);

    end

    varargout{1} = pi*y.^2; % Area
    varargout{2} = 2*pi*y.*dydx; % dAdx
    varargout{3} = 2*y; % D
    varargout{4} = 0.01*ones(length(x),1); % thickness, constant for now
end

end

function y = utoyMap(knots,coefs,i,u)

    y = -(((u - knots(1 + i))*(u - knots(3 + i)))/((knots(1 + i) - knots(3 + i))*(knots(2 + i) - knots(3 + i))) + ((u - knots(2 + i))*(u - knots(4 + i)))/((knots(2 + i) - knots(3 + i))*(knots(2 + i) - knots(4 + i))))*coefs(2, 1 + i) + ((u - knots(3 + i))^2*coefs(2, i))/((knots(1 + i) - knots(3 + i))*(knots(2 + i) - knots(3 + i))) + ((u - knots(2 + i))^2*coefs(2, 2 + i))/((knots(2 + i) - knots(3 + i))*(knots(2 + i) - knots(4 + i)));

end

function x = utoxMap(knots,coefs,i,u)

    x = - (((u - knots(1 + i))*(u - knots(3 + i)))/((knots(1 + i) - knots(3 + i))*(knots(2 + i) - knots(3 + i))) + ((u - knots(2 + i))*(u - knots(4 + i)))/((knots(2 + i) - knots(3 + i))*(knots(2 + i) - knots(4 + i))))*coefs(1, 1 + i) + ((u - knots(3 + i))^2*coefs(1, i))/((knots(1 + i) - knots(3 + i))*(knots(2 + i) - knots(3 + i))) + ((u - knots(2 + i))^2*coefs(1, 2 + i))/((knots(2 + i) - knots(3 + i))*(knots(2 + i) - knots(4 + i)));

end

function result = dgdu(knots,coefs,i,u)

    result = -(2*(u*knots(1 + i)*coefs(1, 1 + i) - u*knots(2 + i)*coefs(1, i) - u*knots(1 + i)*coefs(1, 2 + i) + u*knots(2 + i)*coefs(1, 1 + i) - u*knots(3 + i)*coefs(1, 1 + i) + u*knots(4 + i)*coefs(1, i) + u*knots(3 + i)*coefs(1, 2 + i) - u*knots(4 + i)*coefs(1, 1 + i) - knots(1 + i)*knots(2 + i)*coefs(1, 1 + i) + knots(1 + i)*knots(2 + i)*coefs(1, 2 + i) + knots(2 + i)*knots(3 + i)*coefs(1, i) - knots(2 + i)*knots(3 + i)*coefs(1, 2 + i) - knots(3 + i)*knots(4 + i)*coefs(1, i) + knots(3 + i)*knots(4 + i)*coefs(1, 1 + i)))/((knots(1 + i) - knots(3 + i))*(knots(2 + i) - knots(3 + i))*(knots(2 + i) - knots(4 + i)));

end

function result = g(knots,coefs,i,u,c)

    result = utoxMap(knots,coefs,i,u) - c;

end

function result = dydu(knots,coefs,i,u)

    result = -(2*(u*knots(1 + i)*coefs(2, 1 + i) - u*knots(2 + i)*coefs(2, i) - u*knots(1 + i)*coefs(2, 2 + i) + u*knots(2 + i)*coefs(2, 1 + i) - u*knots(3 + i)*coefs(2, 1 + i) + u*knots(4 + i)*coefs(2, i) + u*knots(3 + i)*coefs(2, 2 + i) - u*knots(4 + i)*coefs(2, 1 + i) - knots(1 + i)*knots(2 + i)*coefs(2, 1 + i) + knots(1 + i)*knots(2 + i)*coefs(2, 2 + i) + knots(2 + i)*knots(3 + i)*coefs(2, i) - knots(2 + i)*knots(3 + i)*coefs(2, 2 + i) - knots(3 + i)*knots(4 + i)*coefs(2, i) + knots(3 + i)*knots(4 + i)*coefs(2, 1 + i)))/((knots(1 + i) - knots(3 + i))*(knots(2 + i) - knots(3 + i))*(knots(2 + i) - knots(4 + i)));

end

function result = dudx(knots,coefs,i,u)

    result = 1/(-(2*(u*knots(1 + i)*coefs(1, 1 + i) - u*knots(2 + i)*coefs(1, i) - u*knots(1 + i)*coefs(1, 2 + i) + u*knots(2 + i)*coefs(1, 1 + i) - u*knots(3 + i)*coefs(1, 1 + i) + u*knots(4 + i)*coefs(1, i) + u*knots(3 + i)*coefs(1, 2 + i) - u*knots(4 + i)*coefs(1, 1 + i) - knots(1 + i)*knots(2 + i)*coefs(1, 1 + i) + knots(1 + i)*knots(2 + i)*coefs(1, 2 + i) + knots(2 + i)*knots(3 + i)*coefs(1, i) - knots(2 + i)*knots(3 + i)*coefs(1, 2 + i) - knots(3 + i)*knots(4 + i)*coefs(1, i) + knots(3 + i)*knots(4 + i)*coefs(1, 1 + i)))/((knots(1 + i) - knots(3 + i))*(knots(2 + i) - knots(3 + i))*(knots(2 + i) - knots(4 + i))));

end

