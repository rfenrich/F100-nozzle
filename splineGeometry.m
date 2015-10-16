function [ varargout ] = splineGeometry( x, outputResult, pp )
% General usage is:
% [output] = nozzleGeometry(x, outputResult, ...) where:
%  x = dimensional location along nozzle length (scalar or vector)
% outputResult = string noting which result is desired to be output (can 
%       take the value of 'A', 'dAdx', 'y', 'D', or 't'. If
%       outputResult is empty, then A, dAdx, D, and t are all output.
%
% For spline parameterized geometry:
% [output] = nozzleGeometry(x, outputResult, pp) where
% pp  = piecewise polynomial coefficient matrix given by the Matlab spline 
% function
% Furthermore, the function will provide the nozzle throat location(s) and
% size(s) when called in the following way:
% [xThroat, yThroat] = nozzleGeometry(x, 'throat', pp)
%
% Rick Fenrich 10/13/15

%y = ppval(pp,x); % cubic interpolant value at x

y = zeros(length(x),1); % initialize
dydx = zeros(length(x),1); % initialize

for ii = 1:length(x)
    ftemp = find(x(ii) <= pp.breaks);
    segment = ftemp(1) - 1; % segment numer of spline
    if segment < 1
        segment = 1;
    end
    y(ii) = pp.coefs(segment,4) + pp.coefs(segment,3)*(x(ii)-pp.breaks(segment)) + pp.coefs(segment,2)*(x(ii)-pp.breaks(segment))^2 + + pp.coefs(segment,1)*(x(ii)-pp.breaks(segment))^3;
    dydx(ii) = pp.coefs(segment,3) + 2*pp.coefs(segment,2)*(x(ii)-pp.breaks(segment)) + 3*pp.coefs(segment,1)*(x(ii)-pp.breaks(segment))^2;
end

D = 2*y;
A = pi*y.^2;
dAdx = 2*pi*y.*dydx;

% For right now, set a constant thickness for the nozzle, regardless of
% shape. Can incorporate variable thickness calculations later above:
t = 0.01*ones(length(x),1);

if(strcmp(outputResult,'A'))
    varargout{1} = A;
elseif(strcmp(outputResult,'dAdx'))
    varargout{1} = dAdx;
elseif(strcmp(outputResult,'y'))
    varargout{1} = y;
elseif(strcmp(outputResult,'D'))
    varargout{1} = D;
elseif(strcmp(outputResult,'t'))
    varargout{1} = t;
elseif(strcmp(outputResult,'throat')) % find location & height of throat for spline only
    
    nn = 0; % counter used for adding values to xThroat, yThroat arrays
    
    % solve quadratic equation to determine where slope = 0 in each cubic
    % polynomial segment of the spline: 1st half of the solutions
    xtemp = (-1/3)*(pp.coefs(:,2)./pp.coefs(:,1)) + sqrt(4*pp.coefs(:,2).^2 - 12*pp.coefs(:,3).*pp.coefs(:,1))./(6*pp.coefs(:,1));
    for ii = 1:length(xtemp)
        if(isreal(xtemp(ii))) % check if root is real
            if(xtemp(ii) + pp.breaks(ii) >= pp.breaks(ii) && xtemp(ii) + pp.breaks(ii) <= pp.breaks(ii+1)) % check if root occurs within spline segment domain
                % calculate value of 2nd derivative at current x location
                ctemp1 = 2*pp.coefs(ii,2) + 6*pp.coefs(ii,1)*xtemp(ii);
                if(ctemp1 > 0) % check if 2nd derivative is positive, implying throat
                    nn = nn + 1;
                    xThroat(nn) = xtemp(ii) + + pp.breaks(ii);
                    yThroat(nn) = pp.coefs(ii,4) + pp.coefs(ii,3).*xtemp(ii) + pp.coefs(ii,2).*xtemp(ii)^2 + pp.coefs(ii,1).*xtemp(ii)^3;
                end
            end
        end
    end

    % solve quadratic equation to determine where slope = 0 in each cubic
    % polynomial segment of the spline: 2nd half of the solutions
    xtemp = (-1/3)*(pp.coefs(:,2)./pp.coefs(:,1)) - sqrt(4*pp.coefs(:,2).^2 - 12*pp.coefs(:,3).*pp.coefs(:,1))./(6*pp.coefs(:,1));
    xtemp = xtemp + pp.breaks(1:end-1)'; % shift x-coordinates for each spline segment
    for ii = 1:length(xtemp)
        if(isreal(xtemp(ii))) % check if root is real
            if(xtemp(ii) + pp.breaks(ii) >= pp.breaks(ii) && xtemp(ii) + pp.breaks(ii) <= pp.breaks(ii+1)) % check if root occurs within spline segment domain
                % calculate value of 2nd derivative at current x location
                ctemp1 = 2*pp.coefs(ii,2) + 6*pp.coefs(ii,1)*xtemp(ii);
                if(ctemp1 > 0) % check if 2nd derivative is positive, implying throat
                    nn = nn + 1;
                    xThroat(nn) = xtemp(ii) + + pp.breaks(ii);
                    yThroat(nn) = pp.coefs(ii,4) + pp.coefs(ii,3).*xtemp(ii) + pp.coefs(ii,2).*xtemp(ii)^2 + pp.coefs(ii,1).*xtemp(ii)^3;
                end
            end
        end
    end
    
    varargout{1} = xThroat;
    varargout{2} = yThroat; % throat diameter
else
    varargout{1} = A;
    varargout{2} = dAdx;
    varargout{3} = D;
    varargout{4} = t;
end

end

