function [ status, shock ] = nozzleState( fluid, pressureRatio, throat, exit )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

gam = fluid.gam;

% Area-Mach function from mass 1-D mass conservation equations:
AreaMachFunc = @(g,M) ((g+1)/2)^((g+1)/(2*(g-1)))*M./(1+(g-1)*M.^2/2).^((g+1)/(2*(g-1)));


% ==================== DETERMINE IDEAL NOZZLE FLOW =======================
options.Display='none'; % used for fsolve
shock.present = false;
shock.M = 0;

exit.critical.Msubsonic = fsolve(@(x) AreaMachFunc(gam,x)*sqrt(throat.Tstag/exit.Tstag)*(exit.Pstag/throat.Pstag) - throat.A/exit.A, 0.5, options);
exit.critical.Msupersonic = fsolve(@(x) AreaMachFunc(gam,x)*sqrt(throat.Tstag/exit.Tstag)*(exit.Pstag/throat.Pstag) - throat.A/exit.A, 2, options);
exit.critical.PtRatioSubsonic = (1 + (gam-1)*exit.critical.Msubsonic^2/2)^(gam/(gam-1));
exit.critical.PtRatioSupersonic = (1 + (gam-1)*exit.critical.Msupersonic^2/2)^(gam/(gam-1));

MbehindShock = sqrt((1 + (gam-1)*exit.critical.Msupersonic^2/2)/(gam*exit.critical.Msupersonic^2 - (gam-1)/2));
exit.normalShock.PtRatio = (exit.A/throat.A)*((gam+1)/2)^((gam+1)/(2*(gam-1)))*MbehindShock*sqrt(1 + (gam-1)*MbehindShock^2/2);

deltaPtRatio = 0.05; % a tolerance for deciding whether fully expanded flow occurs

if (pressureRatio <= 1)
    status = 'no flow';
elseif (pressureRatio < exit.critical.PtRatioSubsonic)
    %fprintf('Subsonic flow throughout\n');
    status = 'subsonic';
    exit.M = sqrt(2/(gam-1))*sqrt(pressureRatio^((gam-1)/gam) - 1);
elseif (pressureRatio < exit.normalShock.PtRatio)
    %fprintf('Shock in nozzle\n');
    status = 'shock';
    exit.M = fsolve( @(x) ((gam+1)/2)^((gam+1)/(2*(gam-1)))*x*sqrt(1 + (gam-1)*x^2/2) - pressureRatio*(throat.A/exit.A)*sqrt(exit.Tstag/throat.Tstag), 0.5, options);
    PtRatio = (throat.A/exit.A)*sqrt(exit.Tstag/throat.Tstag)/AreaMachFunc(gam,exit.M);
    shock.M = fsolve( @(x) (((gam+1)*x^2/2)/(1 + (gam-1)*x^2/2))^(gam/(gam-1))*(((gam+1)/2)/(gam*x^2 - (gam-1)/2))^(1/(gam-1)) - PtRatio,2,options);
    shock.Mpost = sqrt((1 + ((gam-1)/2)*shock.M^2)/(gam*shock.M^2 - (gam-1)/2));
    shock.present = true;
elseif (pressureRatio < exit.critical.PtRatioSupersonic - deltaPtRatio)
    %fprintf('Overexpanded flow\n');
    status = 'overexpanded';
elseif (pressureRatio < exit.critical.PtRatioSupersonic + deltaPtRatio)
    %fprintf('Approximately fully expanded flow\n');
    status = 'fully expanded';
else
    %fprintf('Underexpanded flow\n');
    status = 'underexpanded';
end

end

