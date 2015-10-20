function [varargout] = nozzleGeometry(x, outputResult, Dinlet, nozzleLength, xThroat,areaRatio_inlet_to_throat, areaRatio_outlet_to_throat, shape)
% nozzleGeometry.m is used to provide nozzle radius, diameter, area, change
% in area, and thickness. Use consistent units.
%
% General usage is:
% [output] = nozzleGeometry(x,outputResult,Dinlet,nozzleLength,xThroat,ARatio1,ARatio2,shape)
% where:
%  x = dimensional location along nozzle length (scalar or vector)
% outputResult = string noting which result is desired to be output (can 
%       take the value of 'A', 'dAdx', 'y', 'D', or 't'. If
%       outputResult is empty, then A, dAdx, D, and t are all output.
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
% Rick Fenrich 6/26/15 updated 10/16/15

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
else
    error('shape not supported by nozzleGeometry.m');
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
else
    varargout{1} = A;
    varargout{2} = dAdx;
    varargout{3} = D;
    varargout{4} = t;
end

end

