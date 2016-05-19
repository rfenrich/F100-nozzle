function [ Nf ] = estimateNf( maxTempVec, maxStressVec, deltat )
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
% Note if maxTempVec and maxStressVec are vectors then an Nf value will be
% calculate for each pair of elements from the vectors.
%
% INPUTS:
% maxTempVec = a vector or scalar temperature (K) to use in estimate of Nf
% maxStressVec = a vector or scalar stress (Pa) to use in estimate of Nf
% deltat = the period in hours of each fatigue cycle (required since Nf due
%       to fatigue is based on an S-n diagram, but Nf due to oxidation and
%       creep is calculated from total lifetime hours --> how long is each
%       mechanical fatigue cycle?)
%
% OUTPUTS:
% Nf = a structure containing cycles to failure for fatigue, oxidation,
%       creep, and combined (total)
%
% Rick Fenrich 3/21/16

Nf.total = zeros(length(maxTempVec),1);
Nf.fatigue = zeros(length(maxTempVec),1);
Nf.oxidation = zeros(length(maxTempVec),1);
Nf.creep = zeros(length(maxTempVec),1);

for mm = 1:length(maxTempVec)
    
    maxTemp = maxTempVec(mm);
    maxStress = maxStressVec(mm);

    % Assume minimum stress is 0.
    stressRange = maxStress;
    if(stressRange < 0) % i.e. compression
      stressRange = -stressRange;
      fprintf('sign of stressRange flipped\n');
    end
    if(~isreal(stressRange))
      stressRange = abs(real(stressRange));
      fprintf('\n!!!!!!!!!!!! stressRange is imaginary !!!!!!!!!!!\n');
    end 

    % ===================== ESTIMATE Nf DUE TO FATIGUE =======================
    % Oxidation has a small effect on residual strength; it is neglected here.
    % Fatigue life is estimated (quite poorly) using only two data points and
    % assuming linear behavior on a log-log S-n plot with no lower fatigue
    % limit.

    data.fatigue.stressRange = ([69; 41] - 14)*1e6; % stress range in Pa
    data.fatigue.Nf = [13880; 120000];
    %loglog(data.fatigue.Nf,data.fatigue.deltaSigma);
    Nf.fatigue(mm) = exp(interp1(log(data.fatigue.stressRange),log(data.fatigue.Nf),log(stressRange),'linear','extrap'));

    % ==================== ESTIMATE Nf DUE TO OXIDATION ======================
    % The only data with explicit temperature dependence is for flexural
    % strength, which may be higher than the yield strength of this material.
    % Nevertheless, for a given stress and temperature, the lifetime in hours
    % is extracted from a table using interpolation.

    % The following are B-spline fits to the experimental data. Anything beyond
    % 4000 hours is extrapolated and may be nonsense.
    knots773 = [0 0 0 1 2 3 4 5 6 7 7 7]';
    coefs773 = [0 0 2200 3200 4200 8000 25000 28000 28000; % 26654
             216 216 216 216 216 194.6 73.4 73.4 73.4];
    knots873 = [0 0 0 1 2 3 4 5 6 7 8 8 8]';
    coefs873 = [0 0 500 1000 2000 4000 5000 8000 28000 28000;
             216 216 216 212.8469 120.8990 73.4 73.4 73.4 73.4 73.4];
    knots973 = [0 0 0 1 2 3 4 5 6 7 7 7]';
    coefs973 = [0 0 500 1000 2000 4000 8000 28000 28000;
             216 216 127.9062 79.1310 73.4 73.4 73.4 73.4 73.4];
    % Use polymer rule of thumb to correct flexural strength to yield strength
    coefs773(2,:) = coefs773(2,:)/1.5;
    coefs873(2,:) = coefs873(2,:)/1.5;
    coefs973(2,:) = coefs973(2,:)/1.5;

    % figure; hold on;
    % tPlot = linspace(0,40000,1000)';
    % plot(tPlot,BsplineGeometryMex(tPlot,3,knots773,coefs773')/2);
    % plot(tPlot,BsplineGeometryMex(tPlot,3,knots873,coefs873')/2);
    % plot(tPlot,BsplineGeometryMex(tPlot,3,knots973,coefs973')/2);
    % legend('773','873','973')

    if(maxStress/1e6 >= max(coefs773(2,:)) - 1e-6) % max stress exceeded
        Nf.oxidation(mm) = 0; % i.e. 0
    elseif(maxStress/1e6 < min(coefs773(2,:)) + 1e-6) % assume a fatigue limit
        Nf.oxidation(mm) = 1e9/deltat; % i.e. 1 billion hours or infinity
    elseif(maxTemp < 273.15) 
        Nf.oxidation(mm) = 1e9/deltat; % i.e. 1 billion hours or infinity
    else
        Tm = maxTemp;
        % If temperature range exceeded, fit data to nearest temperature
        if(Tm > 973.15), Tm = 973.15; end

        % Set upper/lower functions for interpolation
        if(Tm <= 773.15)
            fu = @(T) max(coefs773(2,:));
            fl = @(T) BsplineGeometryMex(T,3,knots773,coefs773')/2;
            Tu = 273.15; 
            Tl = 773.15;
            zeroFunc = @(T) maxStress/1e6 - fl(T) - (fu(T) - fl(T))*(Tm - Tl)/(Tu - Tl);
            startGuess = 5000;
        elseif(Tm <= 873.15)
            fu = @(T) BsplineGeometryMex(T,3,knots773,coefs773')/2;
            fl = @(T) BsplineGeometryMex(T,3,knots873,coefs873')/2;
            Tu = 773.15;
            Tl = 873.15;
            zeroFunc = @(T) maxStress/1e6 - fl(T) - (fu(T) - fl(T))*(Tm - Tl)/(Tu - Tl);
            startGuess = 5000;
        elseif(Tm <= 973.15)
            fu = @(T) BsplineGeometryMex(T,3,knots873,coefs873')/2;
            fl = @(T) BsplineGeometryMex(T,3,knots973,coefs973')/2;
            Tu = 873.15;
            Tl = 973.15;
            zeroFunc = @(T) maxStress/1e6 - fl(T) - (fu(T) - fl(T))*(Tm - Tl)/(Tu - Tl);
            % Find sufficient start guess
            Tstart = [5 50 100 1000 5000];
            val = zeros(length(Tstart),1);
            for ii = 1:length(Tstart)
                val(ii) = zeroFunc(Tstart(ii));
            end
            ind = find(val > 0,1,'first');
            startGuess = Tstart(ind);
            %startGuess = 5000;
				else
					fprintf('  ## ERROR in EstimateNf. Return 0.')
					Nf = 0;
					return;
        end

        if(Tm <= 773.15)
            sigmaLimit = min(coefs773(2,:)) + (max(coefs773(2,:)) - min(coefs773(2,:)))*(Tm - 773.15)/(273.15 - 773.15);
            if(maxStress/1e6 < sigmaLimit)
                Nf.oxidation(mm) = 1e9/deltat; % i.e. 1 billion hours
            else
                Nf.oxidation(mm) = fzero(zeroFunc,startGuess)/deltat;
            end
        else

            Nf.oxidation(mm) = fzero(zeroFunc,startGuess)/deltat;
            if(Nf.oxidation(mm) < 0 || Nf.oxidation(mm) > 1e9)
                fprintf('Nf.oxidation < 0\n');
                fprintf('Max stress: %f\n',maxStress);
                fprintf('Max temp: %f\n',maxTemp);
                Nf.oxidation(mm) = 0;
            end
        end

    end

    % ====================== ESTIMATE Nf DUE TO CREEP ========================
    % The data for creep is given only for a temperature of 566 degC. 

    data.creep.stressRoomTemp = [55, 55];
    data.creep.strainRateRoomTemp = [1.5e-9, 1.1e-9];
    data.creep.stress = [69, 83, 96]; % MPa
    data.creep.strainRate = [1.5e-9, 9.3e-10, 3.9e-9]; % s^-1
    data.creep.life = [133, 76, 17.5]; % hr
    %loglog(data.creep.stress,data.creep.strainRate);
    %coefs = polyfit(log(data.creep.stress),log(data.creep.strainRate),1);
    slopeCreep = 2.65; % obtained from linear fit of above data (slope on log-log plot)
    interceptCreep = -31.8; % obtained from " "

    % Assumed value for activation energy for creep
    creepActivationEnergy = 50e3; % J/mol
    coefCreep = exp(interceptCreep + creepActivationEnergy/8.3144598/(566+273.15));

    creepRate = coefCreep*(maxStress/1e6)^slopeCreep*exp(-creepActivationEnergy/8.3144598/maxTemp);
    %fprintf('Creep rate: %e\n',creepRate);
    %if(maxStress/1e6 > 76.17) % estimate creep using linear fit
    %    creepRate = exp(slopeCreep*log(maxStress/1e6) + interceptCreep);
    %end%else
    %    creepRate = 1.5e-9; % s^-1
    %end

    % Chawla CMC textbook assumes failure at 1% strain; experimental data shows
    % failure occuring between 0.1% and 0.5% strain
    creepFailureStrain = 0.003; % assume failure occurs at 0.3%
    T_creep = creepFailureStrain/creepRate/3600; % hr
    Nf.creep(mm) = T_creep/deltat;

    % ========================= ESTIMATE TOTAL Nf ============================
    if(Nf.fatigue(mm) == 0 || Nf.oxidation(mm) == 0 || Nf.creep(mm) == 0)
        Nf.total(mm) = 0;
    else
        Nf.total(mm) = 1/(1/Nf.fatigue(mm) + 1/Nf.oxidation(mm) + 1/Nf.creep(mm));
    end

    % =========================== PRINT RESULTS ==============================
    % fprintf('Nf.fatigue: %e\n',Nf.fatigue);
    % fprintf('Nf.oxidation: %e\n',Nf.oxidation);
    % fprintf('Nf.creep: %e\n',Nf.creep);
    % fprintf('Total Nf: %e\n',Nf.total);
    
end

end

