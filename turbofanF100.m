function [ thrust, sfc, thermalEfficiency, engine ] = turbofanF100( altitude, mach, control, error )
% Calculate performance for F100-PW-220 turbofan engine w/o afterburning.
%
% INPUTS:
% altitude = altitude in ft
% mach = freestream Mach number
% control = control structure dictating which parameters are to externally
%           set by the user; if left empty, standard values are assumed by
%           the program
%
% OUTPUTS:
% thrust = structure with calculated thrust
% sfc = specific fuel consumption
% thermal efficiency of engine
% engine = structure with nozzle parameters
%
% Rick Fenrich 7/31/15 adapted 9/24/15 for error sensitivity study

options.Display='none'; % used for Matlab's fzero

% If control is empty, set all parameters to zero
if(isempty(control))
    control.bypassRatio = 0;
    control.f = 0;
    control.fan.PstagRatio = 0;
    control.fan.efficiency.polytropic = 0;
    control.compressor.efficiency.polytropic = 0;
    control.compressor.overallPressureRatio = 0;
    control.burner.efficiency = 0;
    control.turbine.efficiency.polytropic = 0;
    control.nozzle.inlet.Abypass2Acore = 0;
    control.nozzle.inlet.D = 0;
    control.nozzle.throat.A = 0;
    control.nozzle.Ainlet2Athroat = 0;
    control.nozzle.Aexit2Athroat = 0;
end
    
% ============================== INPUTS ==================================

% Mixing options: 'area-averaged', 'massflow-averaged', 'multistream-samePt'
mixing = 'area-averaged';
%fprintf('%s mixing\n',mixing);

atm = StndAtm(altitude*0.3048,'SI');
freestream.P = atm.P; % Pa, atmospheric pressure
freestream.T = atm.T; % K, atmospheric temperature

gam = 1.4; % ratio of specific heats
Cp = 1006; % J/kg*K, specific heat for dry air
R = 287.06; % J/kg*K, specific gas constant for dry air

freestream.M = mach;
freestream.U = freestream.M*sqrt(gam*R*freestream.T);

% General engine parameters (found in literature)
turbine.inlet.TstagLimit = 1672.15; % K, static temperature at turbine inlet
turbine.inlet.Tstag = turbine.inlet.TstagLimit;
fuel.enthalpy = 4.28e7; % J/kg, Jet A/B, JP-4, or JP-8

% ------------------------ CONTROLLED INPUTS -----------------------------

if(control.bypassRatio)
    bypassRatio = control.bypassRatio;
else
    bypassRatio = 0.6;
end

if(control.f) % fuel mass fraction ratio
    f = control.f; 
else % f will be calculated later
end

% Fan
if(control.fan.PstagRatio)
    fan.PstagRatio = control.fan.PstagRatio;
else
    fan.PstagRatio = 3.06;
end

if(control.fan.efficiency.polytropic)
    fan.efficiency.polytropic = control.fan.efficiency.polytropic;
else
    fan.efficiency.polytropic = 0.83; % includes fan for bypass air
end

% Compressor
if(control.compressor.efficiency.polytropic)
    compressor.efficiency.polytropic = control.compressor.efficiency.polytropic;
else
    compressor.efficiency.polytropic = 0.87; %0.895 % includes compressor for core air
end

if(control.compressor.overallPressureRatio)
    compressor.overallPressureRatio = control.compressor.overallPressureRatio;
else
    compressor.overallPressureRatio = 24.5;
end

% Burner
if(control.burner.PstagRatio)
    burner.PstagRatio = control.burner.PstagRatio;
else
    burner.PstagRatio = 0.95;
end
    
if(control.burner.efficiency)
    burner.efficiency = control.burner.efficiency;
else
    burner.efficiency = 0.95;
end

% Turbine
if(control.turbine.efficiency.polytropic)
    turbine.efficiency.polytropic = control.turbine.efficiency.polytropic;
else
    turbine.efficiency.polytropic = 0.85;
end

if(control.turbine.efficiency.shaft)
    turbine.efficiency.shaft = control.turbine.efficiency.shaft;
else
    turbine.efficiency.shaft = 0.97;
end

% Nozzle geometry

if(control.nozzle.inlet.Abypass2Acore)
    nozzle.inlet.Abypass2Acore = control.nozzle.inlet.Abypass2Acore;
else
    nozzle.inlet.Abypass2Acore = 0.2962;
end

if(control.nozzle.inlet.D)
    nozzle.inlet.D = control.nozzle.inlet.D;
else
    nozzle.inlet.D = 0.651;
end
nozzle.inlet.A = pi*nozzle.inlet.D^2/4;

if(control.nozzle.throat.A)
    nozzle.throat.A = control.nozzle.throat.A;
    if(control.nozzle.exit.A)
        nozzle.exit.A = control.nozzle.exit.A;
    else
        nozzle.exit.A = 0.4207; % exit area is fixed
    end
    nozzle.Ainlet2Athroat = nozzle.inlet.A/nozzle.throat.A;
    nozzle.Aexit2Athroat = nozzle.exit.A/nozzle.throat.A;
elseif(control.nozzle.Ainlet2Athroat && control.nozzle.Aexit2Athroat)
    nozzle.Ainlet2Athroat = control.nozzle.Ainlet2Athroat;
    nozzle.Aexit2Athroat = control.nozzle.Aexit2Athroat;
    nozzle.throat.A = nozzle.inlet.A/nozzle.Ainlet2Athroat;
    nozzle.exit.A = nozzle.Aexit2Athroat*nozzle.inlet.A/nozzle.Ainlet2Athroat;
else
    fprintf('Nozzle area ratios not specified. Setting defaults...\n');
    nozzle.Ainlet2Athroat = 1.368;
    nozzle.Aexit2Athroat = 1.4;
    nozzle.throat.A = nozzle.inlet.A/nozzle.Ainlet2Athroat;
    nozzle.exit.A = nozzle.Aexit2Athroat*nozzle.inlet.A/nozzle.Ainlet2Athroat;    
end

% -------------------------- NOZZLE GEOMETRY -----------------------------
nozzle.length = 1;
nozzle.shape = 'linear';
nozzle.xThroat = 0.33;
nozzle.xExit = nozzle.length;

% ---------------------- GENERAL ENGINE PARAMETERS -----------------------
% Diffuser
diffuser.TstagRatio = 1;

% -------------------------- USEFUL FUNCTIONS ----------------------------

AreaMachFunc = @(g,M) ((g+1)/2)^((g+1)/(2*(g-1)))*M./(1+(g-1)*M.^2/2).^((g+1)/(2*(g-1)));
massFlowRate = @(Pstag,Area,Tstag,Mach) (gam/((gam+1)/2)^((gam+1)/(2*(gam-1))))*Pstag*Area*AreaMachFunc(gam,Mach)/sqrt(gam*R*Tstag);

% =========================== CALCULATIONS ===============================

% ---------------------------- FREESTREAM --------------------------------

freestream.TstagRatio = (1 + (gam-1)*freestream.M^2/2);
freestream.PstagRatio = freestream.TstagRatio^(gam/(gam-1));
freestream.Tstag = freestream.TstagRatio*freestream.T;
freestream.Pstag = freestream.PstagRatio*freestream.P;

% ----------------------------- DIFFUSER ---------------------------------

if(freestream.M < 1)
    diffuser.PstagRatio = 0.97; % modified from original value of 1
elseif(freestream.M < 5)
    diffuser.PstagRatio = 1 - 0.075*(freestream.M - 1)^1.35;
else
    diffuser.PstagRatio = 800/(freestream.M^4 + 935);
end

% ------------------------------- FAN ------------------------------------

fan.inlet.Tstag = freestream.Tstag*diffuser.TstagRatio;
fan.inlet.Pstag = freestream.Pstag*diffuser.PstagRatio;
fan.TstagRatio = fan.PstagRatio^((gam-1)/gam/fan.efficiency.polytropic);
fan.exit.Pstag = fan.inlet.Pstag*fan.PstagRatio;
fan.exit.Tstag = fan.inlet.Tstag*fan.TstagRatio;
fan.efficiency.isentropic = (fan.PstagRatio^((gam-1)/gam) - 1)/(fan.PstagRatio^((gam-1)/gam/fan.efficiency.polytropic) - 1);

% ---------------------------- COMPRESSOR --------------------------------

compressor.inlet.Tstag = fan.inlet.Tstag*fan.TstagRatio;
compressor.inlet.Pstag = fan.inlet.Pstag*fan.PstagRatio;
compressor.PstagRatio = compressor.overallPressureRatio/(freestream.PstagRatio*diffuser.PstagRatio*fan.PstagRatio);
compressor.TstagRatio = (compressor.PstagRatio^((gam-1)/(gam*compressor.efficiency.polytropic)));
compressor.efficiency.isentropic = (compressor.PstagRatio^((gam-1)/gam) - 1)/(compressor.PstagRatio^((gam-1)/gam/compressor.efficiency.polytropic) - 1);

compressor.overallTemperatureRatio = compressor.TstagRatio*fan.TstagRatio;

% ------------------------------ BURNER ----------------------------------

burner.inlet.Tstag = compressor.inlet.Tstag*compressor.TstagRatio;
burner.inlet.Pstag = compressor.inlet.Pstag*compressor.PstagRatio;
burner.TstagRatio = turbine.inlet.Tstag/burner.inlet.Tstag;

% ----------------------------- FUEL FLOW --------------------------------

if(~exist('f','var')) % Calculate f given max turbine inlet Tstag
    tau_lambda = turbine.inlet.Tstag/fan.inlet.Tstag;
    f = (tau_lambda - compressor.overallTemperatureRatio)/((1+bypassRatio)*(burner.efficiency*fuel.enthalpy/(Cp*fan.inlet.Tstag) - tau_lambda));
else % Calculate turbine inlet Tstag given f
    tau_lambda = (f*burner.efficiency*fuel.enthalpy/(Cp*fan.inlet.Tstag) + compressor.overallTemperatureRatio/(1+bypassRatio))/(f + 1/(1+bypassRatio));
    turbine.inlet.Tstag = tau_lambda*fan.inlet.Tstag;
    if(turbine.inlet.Tstag > turbine.inlet.TstagLimit)
        fprintf('! Turbine inlet max temperature reached. Resetting fuel flow.\n');
        tau_lambda = turbine.inlet.TstagLimit/fan.inlet.Tstag;
        burner.TstagRatio = turbine.inlet.TstagLimit/burner.inlet.Tstag;
        f = (tau_lambda - compressor.overallTemperatureRatio)/((1+bypassRatio)*(burner.efficiency*fuel.enthalpy/(Cp*fan.inlet.Tstag) - tau_lambda));
    end
end

% ------------------------------ TURBINE ---------------------------------

turbine.inlet.Pstag = burner.inlet.Pstag*burner.PstagRatio;
turbine.TstagRatio = 1 - (compressor.overallTemperatureRatio - 1 + bypassRatio*(fan.TstagRatio - 1))/(turbine.efficiency.shaft*tau_lambda*(1+f+bypassRatio*f));
turbine.PstagRatio = turbine.TstagRatio^(gam/(gam-1)/turbine.efficiency.polytropic);
turbine.exit.Pstag = turbine.inlet.Pstag*turbine.PstagRatio;
turbine.exit.Tstag = turbine.inlet.Tstag*turbine.TstagRatio;
turbine.efficiency.isentropic = (turbine.PstagRatio^((gam-1)*turbine.efficiency.polytropic/gam) - 1)/(turbine.PstagRatio^((gam-1)/gam) - 1);

% ------------------------------ MIXING ----------------------------------
%tolerance = 1e-10; % tolerance for error in nozzle inlet Mach number
tolerance = error.betweenIterations.inletMach;
errorNozzleInletMach = 1;
iterationLimit = 10;
counter = 0;

while (abs(errorNozzleInletMach) > tolerance)
    
    counter = counter + 1; 

    if (strcmp(mixing,'area-averaged')) 
    % static pressure and temperature area averaged

        if(counter == 1) % first time around must estimate turbine & fan exit Mach numbers
            fan.exit.M = 0.6; % first guess
            incorrect = 1;
            while incorrect % iterate fan.exit.M to satisfy mass conservation equation if necessary
            % Calculate turbine exit Mach number using mass conservation
            [turbine.exit.M, ~, exitflag] = fzero( @(x) AreaMachFunc(gam,fan.exit.M) - bypassRatio*sqrt(fan.exit.Tstag/turbine.exit.Tstag)*(turbine.exit.Pstag/fan.exit.Pstag)*(1/nozzle.inlet.Abypass2Acore)*AreaMachFunc(gam,x),0.5,options);
                if (exitflag ~= 1)
                    %fprintf('! fzero calculated wrong turbine.exit.M, lowering fan.exit.M\n');
                    fan.exit.M = fan.exit.M - 0.1; % lower the fan exit Mach
                else
                    incorrect = 0;
                end
            end
        end

        turbine.exit.P = turbine.exit.Pstag/(1 + (gam-1)*turbine.exit.M^2/2)^(gam/(gam-1));
        turbine.exit.T = turbine.exit.Tstag/(1 + (gam-1)*turbine.exit.M^2/2);

        fan.exit.P = fan.exit.Pstag/(1 + (gam-1)*fan.exit.M^2/2)^(gam/(gam-1));
        fan.exit.T = fan.exit.Tstag/(1 + (gam-1)*fan.exit.M^2/2);

        nozzle.inlet.T = (fan.exit.T)*(nozzle.inlet.Abypass2Acore/(1 + nozzle.inlet.Abypass2Acore)) + (turbine.exit.T)*(1/(1 + nozzle.inlet.Abypass2Acore));
        nozzle.inlet.P = (fan.exit.P)*(nozzle.inlet.Abypass2Acore/(1 + nozzle.inlet.Abypass2Acore)) + (turbine.exit.P)*(1/(1 + nozzle.inlet.Abypass2Acore));    

        if(counter == 1) % first time around, must estimate nozzle inlet Mach number
            nozzle.inlet.M = (fan.exit.M*sqrt(gam*R*fan.exit.T)*(nozzle.inlet.Abypass2Acore/(1 + nozzle.inlet.Abypass2Acore)) + turbine.exit.M*sqrt(gam*R*turbine.exit.T)*(1/(1 + nozzle.inlet.Abypass2Acore)))/sqrt(gam*R*nozzle.inlet.T);
        end

        nozzle.inlet.Tstag = (1 + (gam-1)*nozzle.inlet.M^2/2)*nozzle.inlet.T;
        nozzle.inlet.Pstag = (1 + (gam-1)*nozzle.inlet.M^2/2)^(gam/(gam-1))*nozzle.inlet.P;

    elseif (strcmp(mixing,'massflow-averaged')) 
    % Pt and Tt are averaged by mass flow, ignoring fuel mass flow rate

        error('not yet implemented');
        %nozzle.inlet.Tstag = bypassRatio*fan.exit.Tstag/(1 + bypassRatio) + turbine.exit.Tstag/(1 + bypassRatio);
        %nozzle.inlet.Pstag = bypassRatio*fan.exit.Pstag/(1 + bypassRatio) + turbine.exit.Pstag/(1 + bypassRatio);    

    elseif (strcmp(mixing,'multistream-samePt')) 
    % bypass and core streams do not mix, but both exit through same nozzle; 
    % same static pressure at bypass duct and turbine exits

        error('not yet implemented');

    end

    % ------------------------------ NOZZLE ----------------------------------

    %[ nozzleFlow, nozzle, xPosition ] = nozzleIdeal( struct('gam',gam,'R',R), nozzle.inlet, freestream, nozzle);
    %nozzle.PstagRatio = 0.97;
    [ nozzleFlow, nozzle, xPosition ] = nozzleNonIdeal( struct('gam',gam,'R',R), nozzle.inlet, freestream, nozzle, 400, error);
    errorNozzleInletMach = (nozzleFlow.M(1) - nozzle.inlet.M)/nozzleFlow.M(1);
    fprintf('% Error in nozzle inlet Mach: %f\n',errorNozzleInletMach);

    % Set new nozzle inlet Mach number
    nozzle.inlet.M = nozzleFlow.M(1);
    %options = optimset('TolFun',1e-8,'Display','none');
    options = optimset('TolFun',error.solver.inletMach,'Display','none');
    ftemp = fsolve(@FanTurbineExitMachFunc,[fan.exit.M, turbine.exit.M],options);
    fan.exit.M = ftemp(1);
    turbine.exit.M = ftemp(2);
    
    if(counter == iterationLimit)
        fprintf('Nozzle inlet Mach number not converged.\n');
        break;
    end

end

% Produce warning if the unnatural/unfeasible happens
if(fan.exit.M > 1)
    fprintf('Supersonic flow in bypass fan duct: %f Mach\n',fan.exit.M);
end

nozzle.massFlowRate = massFlowRate(nozzle.inlet.Pstag,nozzle.inlet.A,nozzle.inlet.Tstag,nozzleFlow.M(1));

% Nozzle exit parameters
nozzle.exit.Pstag = nozzle.inlet.Pstag*nozzle.PstagRatio;
nozzle.exit.M = nozzleFlow.M(end);
nozzle.exit.T = nozzleFlow.T(end);
nozzle.exit.U = nozzle.exit.M*sqrt(gam*R*nozzle.exit.T);
nozzle.exit.P = nozzle.exit.Pstag/(1 + (gam-1)*nozzle.exit.M^2/2)^(gam/(gam-1));

% ======================= CALCULATE PERFORMANCE ==========================

% Fuel mass flow rate
fuel.massFlowRate = (f/(1+f))*nozzle.massFlowRate;

% Estimate thrust
thrust.total = nozzle.massFlowRate*(nozzle.exit.U - freestream.U) + fuel.massFlowRate + (nozzle.exit.P - freestream.P)*nozzle.exit.A;
thrust.specific = thrust.total/(nozzle.massFlowRate*(1-f/(1+f)));

% Estimate specific fuel consumption
sfc = fuel.massFlowRate/thrust.total * 2.20462 * 3600 / 0.22481;

% Estimate thermal efficiency
thermalEfficiency = (thrust.total*freestream.U + 0.5*nozzle.massFlowRate*(nozzle.exit.U - freestream.U)^2 - 0.5*fuel.massFlowRate*freestream.U^2)/(fuel.massFlowRate*fuel.enthalpy);

% Output data
engine.freestream = freestream;
engine.diffuser = diffuser;
engine.fan = fan;
engine.compressor = compressor;
engine.burner = burner;
engine.fuel = fuel;
engine.turbine = turbine;
engine.nozzle = nozzle;

% ---------------------------- PRINT DATA --------------------------------

% fprintf('%s\n',nozzle.status);
% fprintf('thrust: %f N\n',thrust.total);
% fprintf('sfc: %0.3f lb/lbf/hr \n',sfc);
% fprintf('eta_th: %0.3f \n',thermalEfficiency);

% fprintf('fan.inlet:\n')
% fan.inlet
% fprintf('fan.exit:\n')
% fan.exit
% fprintf('turbine.exit:\n')
% turbine.exit
% fprintf('nozzle.inlet:\n')
% nozzle.inlet

% fprintf('fan.exit.T: %f\n',fan.exit.T);
% fprintf('turbine.exit.T: %f\n',turbine.exit.T);
% fprintf('nozzle.inlet.Tstag: %f\n',nozzle.inlet.Tstag);
% fprintf('nozzle.inlet.Pstag: %f\n',nozzle.inlet.Pstag);
% fprintf('nozzle.massFlowRate: %f\n',nozzle.massFlowRate);
% fprintf('fuel mass flow rate: %f\n',f);

% ============================= PLOTTING =================================

% % Plot stagnation temp. and pressure through engine
% pStation = [freestream.Pstag; fan.inlet.Pstag; compressor.inlet.Pstag; burner.inlet.Pstag; turbine.inlet.Pstag; nozzle.inlet.Pstag; nozzle.exit.Pstag];
% TStation = [freestream.Tstag; fan.inlet.Tstag; compressor.inlet.Tstag; burner.inlet.Tstag; turbine.inlet.Tstag; nozzle.inlet.Tstag; nozzle.exit.Tstag];
% stationNumber = {'0','2','2.5','3','4','5-7','e'};
% figure
% subplot(1,2,1)
% hold on
% plot(pStation/1e3,'-ro')
% title('Stagnation pressure')
% ylabel('kPa')
% set(gca,'XTickLabel',stationNumber);
% grid on
% subplot(1,2,2)
% hold on
% plot(TStation,'-ro')
% title('Stagnation temperature')
% ylabel('K')
% set(gca,'XTickLabel',stationNumber);
% grid on
%
% formatPlot;

function [ ftemp ] = FanTurbineExitMachFunc( itemp )

    fanExitMach = itemp(1);
    turbineExitMach = itemp(2);
    fanExitT = fan.exit.Tstag/(1 + (gam-1)*fan.exit.M^2/2);
    turbineExitT = turbine.exit.Tstag/(1 + (gam-1)*turbine.exit.M^2/2);
    nozzleInletT = (fanExitT)*(nozzle.inlet.Abypass2Acore/(1 + nozzle.inlet.Abypass2Acore)) + (turbineExitT)*(1/(1 + nozzle.inlet.Abypass2Acore));
    
    ftemp(1) = AreaMachFunc(gam,fanExitMach) - bypassRatio*sqrt(fan.exit.Tstag/turbine.exit.Tstag)*(turbine.exit.Pstag/fan.exit.Pstag)*(1/nozzle.inlet.Abypass2Acore)*AreaMachFunc(gam,turbineExitMach);
    ftemp(2) = (fanExitMach*sqrt(gam*R*fanExitT)*(nozzle.inlet.Abypass2Acore/(1 + nozzle.inlet.Abypass2Acore)) + turbineExitMach*sqrt(gam*R*turbineExitT)*(1/(1 + nozzle.inlet.Abypass2Acore)))/sqrt(gam*R*nozzleInletT) - nozzle.inlet.M;
    
end

end

