function [atm] = StndAtm( varargin )
% StndAtm calculates properties of the 1976 U.S. Standard Atmosphere for
% up to altitudes of 230,000 ft. This function is adapted from Ilan Kroo's
% 1995 JavaScript code found on Stanford's AA 241 webpage. 
%
% StndAtm(H) returns a structure with atmospheric properties corresponding
% to altitude H in English units. H can be a vector or scalar.
%
% StndAtm(H, UNITS) returns a structure with atmospheric properties 
% corresponding to altitude H in the units specified by the string UNITS.
%
% StndAtm(H, V, RL) returns a structure with atmospheric properties 
% corresponding to altitude H, velocity V, and characteristic length RL in 
% English units. H can be a vector or scalar. V must be a vector of the
% same length as H, or it may be a scalar. RL must be a scalar.
%
% StndAtm(H, V, RL, UNITS) returns a structure with atmospheric properties 
% corresponding to altitude H, velocity V, and characteristic length RL in 
% the units specifed by the string UNITS.
%
% The string UNITS can be either 'US' or 'SI'.
%
% Property                  Units (US)         Units (SI)
% ------------------------------------------------------------------------
% Altitude                  ft                 m
% Velocity                  ft/sec             m/s
% Characteristic Length     ft                 m
% Temperature               degR               degK
% Density                   sl/ft^3            kg/m^3
% Pressure                  lb/ft^2            N/m^2
% Kinematic Viscosity       ft/sec             m/s
% Dynamic Viscosity         lb-sec/ft^2        N-s/m^2
% Dynamic Pressure          lb/ft^2            N/m^2
% Speed of Sound            ft/sec             m/s
% Mach Number               -                  -
% Reynolds Number           -                  -
% Turbulent Friction Coeff. -                  -
% Laminar Friction Coeff.   -                  -
% Critical Pressure Coeff.  -                  -
% Minimum Pressure Coeff.   -                  -
%
% Rick Fenrich 9/12/14

% ======================= INITIALIZE VARIABLES ===========================

if(nargin == 1)
    atm.h = varargin{1};
    units = 'US';
elseif(nargin == 2)
    atm.h = varargin{1};
    units = varargin{2};
elseif(nargin == 3)
    atm.h = varargin{1};
    atm.V = varargin{2};
    atm.rl = varargin{3};
    units = 'US';
elseif(nargin == 4)
    atm.h = varargin{1};
    atm.V = varargin{2};
    atm.rl = varargin{3};
    units = varargin{4};
elseif(nargin == 5)
    fprintf('StndAtm.m: Ignoring extra input arguments\n');
    atm.h = varargin{1};
    atm.V = varargin{2};
    atm.rl = varargin{3};
    units = varargin{4};
end

h = atm.h;
StndAtmNargin = nargin;

if(StndAtmNargin > 2)
    if(length(atm.V) ~= length(h))
        atm.V = atm.V*ones(length(h),1);
    end
end

atm.T = zeros(length(h),1); % temperature
atm.rho =  zeros(length(h),1); % density
atm.P  = zeros(length(h),1); % pressure
atm.c  = zeros(length(h),1); % speed of sound
atm.mu =  zeros(length(h),1); % dynamic viscosity

if(StndAtmNargin > 2)
    atm.M = zeros(length(h),1); % Mach number
    atm.q  = zeros(length(h),1); % dynamic pressure
    atm.cpstar  = zeros(length(h),1); % critical CP
    atm.cpmin  = zeros(length(h),1); % minimum CP
    atm.Re =  zeros(length(h),1); % Reynolds number
    atm.cflam  = zeros(length(h),1); % laminar coefficient of friction
    atm.cfturb  =  zeros(length(h),1); % turbulent coefficient of friction
end

% ================== CALCULATE ATMOSPHERE PROPERTIES =====================

if(strncmp(units,'US',2)) % US units are given
    
    % Calculate atmospheric properties
    for ii = 1:length(h)
        [atm] = Compute(atm);
    end
    
    atm.units = units;
    
elseif(strncmp(units,'SI',2)) % SI units are given
    
	h = h*3.28084; % convert to ft
    if(StndAtmNargin > 2)
        atm.V = atm.V*3.28084;
        atm.rl = atm.rl*3.28084;
    end
    
    % Calculate atmospheric properties
    for ii = 1:length(h)
        [atm] = Compute(atm);
    end

    % convert to SI
    atm.T = atm.T/1.8;
    atm.rho = atm.rho / .068521  / .028317;
    atm.P = atm.P / .020885;
    atm.c = atm.c/3.2808;
    atm.mu = atm.mu/.22481/.092903;
    if(StndAtmNargin > 2)
        atm.q = atm.q / .020885;
    end
    
    atm.units = units;
    
end

% ======================= ADDITIONAL FUNCTIONS ===========================

function [atm] = Compute(atm)
   
   TEMPSL = 518.67; % degR
   RHOSL = 0.00237689; % slug/ft^3
   PRESSSL = 2116.22; % lbf/ft^2
   saTheta = 1.0;
   saSigma = 1.0;
   saDelta = 1.0;

   if ( h(ii)>232940 )
       error('Atmospheric model only valid up to 232,2940 ft');
   end
   if ( h(ii)<232940 )
      saTheta = 1.434843 - h(ii)/337634;
      saSigma = (0.79899-h(ii)/606330)^11.20114;
      saDelta = (0.838263-h(ii)/577922)^12.20114;
   end
   if ( h(ii)<167323 )
      saTheta = 0.939268;
      saSigma = 0.00116533 * exp( (h(ii)-154200)/-25992 );
      saDelta = 0.00109456 * exp( (h(ii)-154200)/-25992 );
   end
   if ( h(ii)<154199 )
      saTheta = 0.482561 + h(ii)/337634;
      saSigma = (0.857003+h(ii)/190115)^-13.20114;
      saDelta = (0.898309+h(ii)/181373)^-12.20114;
   end
   if ( h(ii)<104987 )
      saTheta = 0.682457 + h(ii)/945374;
      saSigma = (0.978261+h(ii)/659515)^-35.16319;
      saDelta = (0.988626+h(ii)/652600)^-34.16319;
   end
   if ( h(ii)<65617 )
      saTheta = 0.751865;
      saSigma = 0.297076 * exp( (36089-h(ii))/20806 );
      saDelta = 0.223361 * exp( (36089-h(ii))/20806 );
   end
   if ( h(ii)<36089 )
      saTheta = 1.0 - h(ii)/145442;
      saSigma = (1.0-h(ii)/145442)^4.255876;
      saDelta = (1.0-h(ii)/145442)^5.255876;
   end

   atm.T(ii) = TEMPSL * saTheta;
   atm.rho(ii) = RHOSL * saSigma;
   atm.P(ii) = PRESSSL * saDelta;
   atm.mu(ii) = 0.0226968*atm.T(ii)^1.5 / ((atm.T(ii))+198.72) / 1000000.0;
   atm.c(ii) = sqrt( 1.4*1716.56*(atm.T(ii)) );
   
   if(StndAtmNargin > 2)
       atm.M(ii) = atm.V(ii)/atm.c(ii);
       atm.q(ii) = 0.7*atm.P(ii)*atm.M(ii)*atm.M(ii);
       atm.Re(ii) = atm.V(ii)*atm.rl*atm.rho(ii)/atm.mu(ii);
       atm.cfturb(ii) = 0.455/(log(atm.Re(ii))/log(10)^2.58);
       atm.cflam(ii) = 1.328/sqrt(atm.Re(ii));

       atm.cpstar(ii) =  (((1/1.2 + atm.M(ii)*atm.M(ii)/6.0)^3.5)-1)/(0.7*atm.M(ii)*atm.M(ii));
       atm.cpmin(ii) =  -1.0/(0.7*atm.M(ii)*atm.M(ii));
   end

end

end