function [varargout] = nozzleGeometrySpline(x,xNode,yNode,slopes,outputResult)
% Takes parameters for a nozzle geometry and a distance x along the length
% of the nozzle and provides area A, rate of change of A (dA/dx), radius y, 
% diameter D, and wall thickness t. Provide consistent units.
%
% INPUTS:
% x = dimensional location along nozzle length (scalar or vector)
% Dinlet = inlet diameter (scalar)
% nozzleLength = nozzle total length (scalar)
% xThroat = location of throat (scalar)
% ARatio1 = ratio of inlet area to throat area
% ARatio2 = ratio of outlet area to throat area
% shape = string determining shape interpolation from inlet to throat to
% outlet
% outputResult = string noting which result is desired to be output (can 
%       take the value of 'A', 'dAdx', 'y', 'D', or 't'. If
%       outputResult is empty, then A, dAdx, D, and t are all output.
%
% OUTPUTS:
% A = area corresponding to nozzle cross-sections at location x
% dAdx = rate of change of area corresponding to " "
% y = radius corresponding to " "
% D = diameter corresponding to " "
%
% Rick Fenrich 6/26/15

% ================= CALCULATE GEOMETRY BASED ON SHAPE ====================

pp = spline(xNode,[slopes(1); yNode; slopes(2)]); % perform piecewise cubic spline interpolation
y = ppval(pp,x); % cubic interpolant value at x
D = 2*y;
A = pi*y.^2;
for ii = 1:length(x)
    ftemp = find(x(ii) <= pp.breaks);
    segment = ftemp(1) - 1; % segment numer of spline
    if segment < 1
        segment = 1;
    end
    dAdx(ii) = pp.coefs(segment,2) + 2*pp.coefs(segment,3)*x(ii) + 3*pp.coefs(segment,4)*x(ii)^2;
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

