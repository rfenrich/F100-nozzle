function [ volume ] = wallVolume( varargin )
% wallVolume calculates the volume of the wall of an axisymmetric
% hollow shape. It accepts three types of geometry:
% 1) inner wall parameterized by cubic spline, thickness parameterized by
% piecewise linear function
% 2) inner wall parameterized by B-spline, thickness parameterized by
% piecewise linear function
% 3) inner wall and thickness parameterized by arbitrary function, volume
% calculated by very fine trapezoidal integration
%
% INPUTS:
% For a cubic spline, vargin = ( pp, pl ) where:
% pp  = piecewise cubic spline struct given by the Matlab spline 
%       function
% pl = piecewise linear struct
%
% For a B-spline, vargin = ( knots, coefs, wall, 'bSpline' ) where:
% knots = knots vector of B-spline
% coefs = coefficients matrix of B-spline
% wall = structure with 3 fields: shape, seed, and breaks determining the
% shape piece-wise linear thickness function
%
% For calculation of geometry via integration for any geometry
% parameterization, vargin = ( length, D, t, 'integrate' ) where:
% length = nozzle.geometry.length
% D = function D(x) which returns diameter as a function of axial distance
% t = function t(x) which returns thickness as a function of axial dist.
%
% OUTPUTS:
% volume = volume of wall (not volume contained by inner wall of nozzle)
%
% Rick Fenrich modified 5/23/16

if(nargin == 2) % assume inputs correspond to cubic spline
    
    pp = varargin{1};
    pl = varargin{2};
    
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
    
elseif(nargin == 4 && strcmp(varargin{4},'B-spline'))
    
    knots = varargin{1};
    coefs = varargin{2};
    wall = varargin{3};
    
    t = @(x) piecewiseLinearGeometry(x,'t',wall); % m, thickness of wall
    
    nozzle.geometry.shape = 'B-spline-mex';
    nozzle.geometry.bSpline.knots = knots;
    nozzle.geometry.bSpline.coefs = coefs;
    
    % Determine nozzle throat
    nozzle.geometry.bSpline.degree = length(nozzle.geometry.bSpline.knots) - length(nozzle.geometry.bSpline.coefs) - 1;
    if(nozzle.geometry.bSpline.degree == 2)
        [xThroat, yThroat] = BsplineGeometry(0, 'throat', nozzle.geometry.bSpline.knots, nozzle.geometry.bSpline.coefs);
    elseif(nozzle.geometry.bSpline.degree == 3)
        [xThroat, yThroat] = BsplineGeometry3(0, 'throat', nozzle.geometry.bSpline.knots, nozzle.geometry.bSpline.coefs);        
    end
    nozzle.geometry.xThroat = xThroat;

    % Define other geometry parameters
    nozzle.geometry.length = nozzle.geometry.bSpline.coefs(1,end);
    nozzle.inlet.D = nozzle.geometry.bSpline.coefs(2,1)*2;
    nozzle.geometry.Ainlet2Athroat = (nozzle.inlet.D/2)^2/yThroat^2;
    nozzle.geometry.Aexit2Athroat = (nozzle.geometry.bSpline.coefs(2,end))^2/yThroat^2;

    [ ~, ~, D, ~ ] = nozzleParameterization( nozzle );
    
    xVolume = linspace(0,nozzle.geometry.length,500)';
    volumeIntegrand = pi*D(xVolume).*t(xVolume) + pi*t(xVolume).^2;
    volume = (xVolume(2)-xVolume(1))*trapz(volumeIntegrand);    
    
elseif(nargin == 4 && strcmp(varargin{4},'integrate'))
    
    nozzleLength = varargin{1};
    D = varargin{2};
    t = varargin{3};
    
    xVolume = linspace(0,nozzleLength,500)';
    volumeIntegrand = pi*D(xVolume).*t(xVolume) + pi*t(xVolume).^2;
    volume = (xVolume(2)-xVolume(1))*trapz(volumeIntegrand);
    
else
    
    error('these inputs not supported for calculating wall volume');
    
end    

end

