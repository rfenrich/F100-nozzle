function [ volume ] = wallVolume( pp, pl )
% wallVolume calculates the volume of the wall of an axisymmetric
% hollow shape. It assumes the interior of the wall is parameterized by
% piecewise cubic splines, but that the thickness of the wall is
% parameterized by a piecewise linear function. The control points need not
% for the wall interior and the wall thickness need not be the same.
%
% INPUTS:
% pp  = piecewise cubic spline struct given by the Matlab spline 
%       function
% pl = piecewise linear struct
%
% OUTPUTS:
% volume = analytically-derived, exact volume of wall

volume = 0;

rControl = pp.breaks;
tControl = pl.breaks;

% Build segment vectors
rSegment = [1];
tSegment = [1];
xSegment = rControl(1);
rr = 1; % counter for index of radius
tt = 1; % counter for index of thickness
nn = 1; % counter for number of segments
while rr <= length(rControl)-1 && tt <= length(tControl)-1
    
    nn = nn + 1; % next segment
    
   if(rControl(rr+1) > tControl(tt+1))
       tt = tt + 1;
       rSegment(nn) = rr;
       tSegment(nn) = tt;
       xSegment(nn) = tControl(tt);
   elseif(rControl(rr+1) < tControl(tt+1))
       rr = rr + 1;
       rSegment(nn) = rr;
       tSegment(nn) = tt;
       xSegment(nn) = rControl(rr);
   else % equal
       tt = tt + 1;
       rr = rr + 1;       
       rSegment(nn) = rr;
       tSegment(nn) = tt;
       xSegment(nn) = rControl(rr);
   end
    
end
rSegment = rSegment(1:end-1);
tSegment = tSegment(1:end-1);

for ii = 1:length(rSegment)
    
    rs = rSegment(ii); % spline segment
    ts = tSegment(ii); % piecewise linear segment
   
    % Parameters for spline segment
    a = pp.coefs(rs,1);
    b = pp.coefs(rs,2);
    c = pp.coefs(rs,3);
    d = pp.coefs(rs,4);
    x1 = pp.breaks(rs);
    
    % Parameters for piecewise linear segment
    % Check if yControl coordinates need to be interpolated from a seed
    if isequal(pl.breaks, pl.seed(:,1))
        yControl = pl.seed(:,2);
    else
        yControl = interpLinear(pl.seed(:,1),pl.seed(:,2),xControl);
    end
    m = (yControl(ts+1)-yControl(ts))/(tControl(ts+1)-tControl(ts));
    x2 = tControl(ts);
    y2 = yControl(ts);
    
    % X coordinates at beginning and end of segment
    xi = xSegment(ii);
    xo = xSegment(ii+1);
    
    vol_contr = pi*xo*((y2 - m*x2)*(- 2*a*x1^3 + 2*b*x1^2 - 2*c*x1 ...
        + 2*d) + (y2 - m*x2)^2) - pi*xi*((y2 - m*x2)*(- 2*a*x1^3 + ...
        2*b*x1^2 - 2*c*x1 + 2*d) + (y2 - m*x2)^2) - ...
        (pi*xi^3*(m*(6*a*x1^2 - 4*b*x1 + 2*c) + m^2 + ...
        (y2 - m*x2)*(2*b - 6*a*x1)))/3 + (pi*xo^3*(m*(6*a*x1^2 - ...
        4*b*x1 + 2*c) + m^2 + (y2 - m*x2)*(2*b - 6*a*x1)))/3 - ...
        (pi*xi^2*(m*(- 2*a*x1^3 + 2*b*x1^2 - 2*c*x1 + 2*d) + ...
        (y2 - m*x2)*(6*a*x1^2 - 4*b*x1 + 2*c) + 2*m*(y2 - m*x2)))/2 ...
        + (pi*xo^2*(m*(- 2*a*x1^3 + 2*b*x1^2 - 2*c*x1 + 2*d) + ...
        (y2 - m*x2)*(6*a*x1^2 - 4*b*x1 + 2*c) + 2*m*(y2 - m*x2)))/2 ...
        - (pi*xi^4*(m*(2*b - 6*a*x1) + 2*a*(y2 - m*x2)))/4 + ...
        (pi*xo^4*(m*(2*b - 6*a*x1) + 2*a*(y2 - m*x2)))/4 - ...
        (2*pi*a*m*xi^5)/5 + (2*pi*a*m*xo^5)/5;

    volume = volume + vol_contr;
    
end

end

