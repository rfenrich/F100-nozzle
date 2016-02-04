function [ varargout ] = piecewiseLinearGeometry( x, outputResult, pl )
% piecewiseLinearGeometry.m provides the function value t(x) of a piecewise
% linear function at (a) given x location(s). It interpolates to find t(x)
% and the slope dtdx(x). The piecewise linear function is built from a seed
% containing a set of coordinates. The piecewise geometry need not use the 
% seed's coordinates as control points, and alternate control point
% x-positions can be defined.
%
% For example, t(x) can be considered as the thickness of an axisymmetric
% shell's wall. Thus dtdx(x) is the change in thickness along the axis of 
% the shell. 
%
% General usage is:
% [output] = nozzleGeometry(x, outputResult, pl) where:
%  x = dimensional location along nozzle length (scalar or vector)
% outputResult = string noting which result is desired to be output (can 
%       take the value of 't', 'dtdx', 'A', 'dAdx', 'y', or 'D'. If
%       outputResult is empty, then t and dtdx are both output.
% pl = struct containing the following fields:
%      seed = nx2 matrix containing x coordinates and vertical distances of 
%             n control points from reference (i.e. vertically measured 
%             thickness of hollow cylinder)
%      breaks = vector of length n containing desired x 
%             coordinates of control points
%
% OUTPUTS:
% t = thickness corresponding to nozzle cross-sections at location x
% dtdx = rate of change of thickness corresponding to " "
% V = total volume of shell of thickness t(x) with base(x)
%
% Rick Fenrich 10/23/15

% ======================= CALCULATE GEOMETRY =============================

% Initialize arrays
t = zeros(length(x),1);
dtdx = zeros(length(x),1);

xControl = pl.breaks;

% Check if yControl coordinates need to be interpolated from a seed
if isequal(pl.breaks, pl.seed(:,1))
    yControl = pl.seed(:,2);
else
    yControl = interpLinear(pl.seed(:,1),pl.seed(:,2),xControl);
end

for ii = 1:length(x)
    ftemp = find(x(ii) <= xControl);
    if(numel(ftemp) == 0) % out of bounds
        ftemp = xControl(1);
    end
    segment = ftemp(1) - 1; % segment numer of piece
    if segment < 1
        segment = 1;
    end
    b = (yControl(segment+1)-yControl(segment))/(xControl(segment+1)-xControl(segment));
    t(ii) = b*(x(ii)-xControl(segment)) + yControl(segment);
    dtdx(ii) = b;
end

% ========================== OUTPUT RESULTS ==============================

if(strcmp(outputResult,'t'))
    varargout{1} = t;
elseif(strcmp(outputResult,'dtdx'))
    varargout{1} = dtdx;
elseif(strcmp(outputResult,'A'))
    error('not yet implemented');%varargout{1} = A;
elseif(strcmp(outputResult,'dAdx'))
    error('not yet implemented');%varargout{1} = dAdx;
elseif(strcmp(outputResult,'y'))
    error('not yet implemented');%varargout{1} = y;
elseif(strcmp(outputResult,'D'))
    error('not yet implemented');%varargout{1} = D;
else
    varargout{1} = t;
    varargout{2} = dtdx;
end

end

