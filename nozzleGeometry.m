function [varargout] = nozzleGeometry(x, outputResult, varargin)
% nozzleGeometry.m is used to provide nozzle radius, diameter, area, change
% in area, and thickness. Use consistent units. When geometry is
% parameterized by a spline, it can also provide the location and size of
% the nozzle throat(s).
%
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
% For other nozzle geometries:
% [output] = nozzleGeometry(x,outputResult,Dinlet,nozzleLength,xThroat,ARatio1,ARatio2,shape)
% where:
% Dinlet = inlet diameter (scalar)
% nozzleLength = nozzle total length (scalar)
% xThroat = location of throat (scalar)
% ARatio1 = ratio of inlet area to throat area
% ARatio2 = ratio of outlet area to throat area
% shape = string determining shape interpolation from inlet to throat to
% outlet
%
% OUTPUTS:
% A = area corresponding to nozzle cross-sections at location x
% dAdx = rate of change of area corresponding to " "
% y = radius corresponding to " "
% D = diameter corresponding to " "
%
% Rick Fenrich 6/26/15 updated 9/1/15

if(length(varargin) == 1)
    pp = varargin{1};
    shape = 'spline';
elseif(length(varargin) == 6)
    Dinlet = varargin{1};
    nozzleLength = varargin{2};
    xThroat = varargin{3};
    areaRatio_inlet_to_throat = varargin{4};
    areaRatio_outlet_to_throat = varargin{5};
    shape = varargin{6};
else
    error('Incorrect number of arguments handed to nozzleGeometry.');
end

% ================= CALCULATE GEOMETRY BASED ON SHAPE ====================

if(strcmp(shape,'tube')) % for straight tube
    A = pi*Dinlet^2/4;
    dAdx = 0;
    y = Dinlet/2;
elseif(strcmp(shape,'linear')) % For linear shape
    Dthroat = sqrt(Dinlet^2/areaRatio_inlet_to_throat);
    Doutlet = sqrt(Dthroat^2*areaRatio_outlet_to_throat);
    y = zeros(length(x),1);
    A = zeros(length(x),1);
    dAdx = zeros(length(x),1);
    D = zeros(length(x),1);
    for ii = 1:length(x)
        if(x(ii) < xThroat) % i.e. before throat
            m = 0.5*(Dthroat - Dinlet)/xThroat;
            y(ii) = m*(x(ii) - xThroat) + Dthroat/2;
            A(ii) = pi*y(ii)^2;
            dAdx(ii) = 2*pi*y(ii)*m;
        else % i.e. after throat
            m = 0.5*(Doutlet - Dthroat)/(nozzleLength - xThroat);
            y(ii) = m*(x(ii) - nozzleLength) + Doutlet/2;
            A(ii) = pi*y(ii)^2;
            dAdx(ii) = 2*pi*y(ii)*m;    
        end
        D(ii) = 2*y(ii);
    end
elseif(strcmp(shape,'cosine')) % for cosine shape
    amplitude = (Dinlet - sqrt(Dinlet^2/ARatio1))/2;
    A = (pi/4)*(Dinlet - amplitude + amplitude*cos(2*pi*x/nozzleLength)).^2;
    dAdx = -pi^2*amplitude*sin(2*pi*x/nozzleLength).*(Dinlet-amplitude+amplitude*cos(2*pi*x/nozzleLength))/nozzleLength;
    D = sqrt(4*A/pi);
    y = D/2;
elseif(strcmp(shape,'spline')) % for spline
    y = ppval(pp,x); % cubic interpolant value at x
    D = 2*y;
    A = pi*y.^2;
    for ii = 1:length(x)
        ftemp = find(x(ii) <= pp.breaks);
        segment = ftemp(1) - 1; % segment numer of spline
        if segment < 1
            segment = 1;
        end
        dydx(ii) = pp.coefs(segment,3) + 2*pp.coefs(segment,2)*(x(ii)-pp.breaks(segment)) + 3*pp.coefs(segment,1)*(x(ii)-pp.breaks(segment))^2;
    end
    dAdx = 2*pi*y.*dydx;
end

% For right now, set a constant thickness for the nozzle, regardless of
% shape. Can incorporate variable thickness calculations later above:
t = 0.01*ones(length(x),1);

% ========================== OUTPUT RESULTS ==============================

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
    varargout{2} = yThroat*2; % throat diameter
else
    varargout{1} = A;
    varargout{2} = dAdx;
    varargout{3} = D;
    varargout{4} = t;
end

end

