function [ Nf ] = estimateNf( maxTemp, maxStress, deltat )
% estimateNf.m very roughly estimates the number of cycles to failure Nf
% for a 2-D test piece of the SiOC Next 312 ceramic matrix composite
% material. All material data is extracted from the characterization of the
% aforementioned material as found in Bansal's Handbook of Ceramic
% Composites. Data for the 8X infiltration material is used. It cannot be
% over-emphasized how approximate the below calculations are. Although
% there is some use of basic material physics equations, much of the data
% is interpolated, extrapolated to unknown regions and based on very few
% data points. Hopefully, it at least captures some ballpark Nf numbers.
%
% Cycles to failure Nf is calculated as a composite of the cycles to
% failure due to mechanical fatigue (ignoring temperature/oxidation
% effects), due to oxidation, and due to creep (likely the dominant
% mechanism for high stress, high temperature loading).
%
% INPUTS:
% maxTemp = a temperature (K) to use in estimate of Nf (usually the max)
% maxStress = a stress (Pa) to use in estimate of Nf (usually the max, 
%       should probably be colocated with the maxTemp)
% deltat = the period in hours of each fatigue cycle (required since Nf due
%       to fatigue is based on an S-n diagram, but Nf due to oxidation and
%       creep is calculated from total lifetime hours --> how long is each
%       mechanical fatigue cycle?)
%
% OUTPUTS:
% Nf = a structure containing cycles to failure for fatigue, oxidation,
%       creep, and combined (total)
%
% Rick Fenrich 3/16/16

% Assume minimum stress is 0.
stressRange = maxStress;

% ===================== ESTIMATE Nf DUE TO FATIGUE =======================
% Oxidation has a small effect on residual strength; it is neglected here.
% Fatigue life is estimated (quite poorly) using only two data points and
% assuming linear behavior on a log-log S-n plot with no lower fatigue
% limit.

data.fatigue.stressRange = ([69; 41] - 14)*1e6; % stress range in Pa
data.fatigue.Nf = [13880; 120000];
%loglog(data.fatigue.Nf,data.fatigue.deltaSigma);
Nf.fatigue = exp(interp1(log(data.fatigue.stressRange),log(data.fatigue.Nf),log(stressRange),'linear','extrap'));

% ==================== ESTIMATE Nf DUE TO OXIDATION ======================
% The only data with explicit temperature dependence is for flexural
% strength, which may be higher than the yield strength of this material.
% Nevertheless, for a given stress and temperature, the lifetime in hours
% is extracted from a table using interpolation.

data.oxidation.stress = [215.7, 223.5, 244.2, 275.4, 270.4, 226.7, 216.0;
                         215.7, 218.5, 240.1, 229.8, 205.5, 132.6, 75.9;
                         215.7, 222.3, 234.0, 132.8, 87.1, 74.0, 73.4]'; % MPa
data.oxidation.temp = [500, 600, 700]' + 273.15; % K
data.oxidation.time = [0, 24, 100, 500, 1000, 2000, 4000]'; % hours
% figure; hold on;
% for ii = 1:length(data.oxidation.temp)
%     plot(data.oxidation.time,data.oxidation.stress(:,ii));
% end
% legend('773 K','873 K','973 K');

% Modify data.oxidation.stress so that the lifetime of a part that would
% have failed at time 0 does not have a lifetime greater than 0 (due to the
% "hill" in the beginning of the data).
data.oxidation.stress(:,1) = 215.7*ones(length(data.oxidation.stress),1);
data.oxidation.stress(2:4,2) = data.oxidation.stress(1,2);
data.oxidation.stress(2:3,3) = data.oxidation.stress(1,3);
% Modify data.oxidation.stress so effect of oxidation stabilizes
data.oxidation.stress(end+1,:) = data.oxidation.stress(end,:);
data.oxidation.time(end+1) = 8000;

% Interpolate stress-time curve for given temperature
if (maxTemp < min(data.oxidation.temp)) % temp is lower than data gives
    % assume oxidation is negligible, arbitrarily choose:
    Nf.oxidation = 120000;
elseif (maxTemp > max(data.oxidation.temp)) % temp is higher than data gives
    % assume oxidation occurs at rate given by data for 973.15 K
    oxidationStress = data.oxidation.stress(:,3);
    if(maxStress/1e6 < min(oxidationStress)) % assume failure will not occur
        Nf.oxidation = 120000;
    else
        [~,ia,~] = unique(oxidationStress);
        T_oxidation = interp1(oxidationStress(ia),data.oxidation.time(ia),maxStress/1e6,'linear','extrap');
        Nf.oxidation = T_oxidation/deltat;
    end
else % interpolate from data
    oxidationStress = interp1(data.oxidation.temp,data.oxidation.stress',maxTemp,'linear','extrap');
    if(maxStress/1e6 < min(oxidationStress)) % assume failure will not occur
        Nf.oxidation = 120000;
    else
        [~,ia,~] = unique(oxidationStress);
        T_oxidation = interp1(oxidationStress(ia),data.oxidation.time(ia),maxStress/1e6,'linear','extrap');
        Nf.oxidation = T_oxidation/deltat;
    end
    %plot(data.oxidation.time,oxidationStress);
end

% ====================== ESTIMATE Nf DUE TO CREEP ========================
% The data for creep is given only for a temperature of 566 degC. 

data.creep.stress = [69, 83, 96]; % MPa
data.creep.strainRate = [1.5e-9, 9.3e-10, 3.9e-9]; % s^-1
data.creep.life = [133, 76, 17.5]; % hr
%loglog(data.creep.stress,data.creep.strainRate);
%coefs = polyfit(log(data.creep.stress),log(data.creep.strainRate),1);
slopeCreep = 2.65; % obtained from linear fit of above data (slope on log-log plot)
interceptCreep = -31.8; % obtained from " "

if(maxStress/1e6 > 76.17) % estimate creep using linear fit
    creepRate = exp(slopeCreep*log(maxStress/1e6) + interceptCreep);
else
    creepRate = 1.5e-9; % s^-1
end

creepFailureStrain = 0.01; % assume failure occurs at 1% (Chawla CMC textbook)
T_creep = creepFailureStrain/creepRate/3600; % hr
Nf.creep = T_creep/deltat;

% ========================= ESTIMATE TOTAL Nf ============================
Nf.total = 1/(1/Nf.fatigue + 1/Nf.oxidation + 1/Nf.creep);

% =========================== PRINT RESULTS ==============================
% fprintf('Nf.fatigue: %e\n',Nf.fatigue);
% fprintf('Nf.oxidation: %e\n',Nf.oxidation);
% fprintf('Nf.creep: %e\n',Nf.creep);
% fprintf('Total Nf: %e\n',Nf.total);



