clc; clear all; close all;

% Convert tabular data from a centered parameter study into an overlaid set
% of 1D plots.

A = importdata('TurbofanUQ-centered-param.tab');

[val ind] = sort(A(:,1));
A = A(ind,:);

center = A(1,:);

num_vars = 10; num_steps = 1000; num_1d_pts = 2*num_steps+1;

Vars         = zeros(num_1d_pts,num_vars);
thrust       = zeros(num_1d_pts,num_vars);
sfc = zeros(num_1d_pts,num_vars);
massFlow = zeros(num_1d_pts,num_vars);
thermalEff = zeros(num_1d_pts,num_vars);

PS_row = 2;
for var=1:num_vars
    for resp_row=1:num_steps
        Vars(resp_row,var)         = A(PS_row,var+1);
        thrust(resp_row,var)       = A(PS_row,num_vars+2);
        sfc(resp_row,var) = A(PS_row,num_vars+3);
        massFlow(resp_row,var) = A(PS_row,num_vars+4);
        thermalEff(resp_row,var) = A(PS_row,num_vars+5);
        PS_row = PS_row+1;
    end
    resp_row = num_steps+1;
    Vars(resp_row,var)           = center(var+1);
    thrust(resp_row,var)         = center(num_vars+2);
    sfc(resp_row,var)   = center(num_vars+3);
    massFlow(resp_row,var)   = center(num_vars+4);
    thermalEff(resp_row,var)   = center(num_vars+5);
    for resp_row=num_steps+2:2*num_steps+1
        Vars(resp_row,var)         = A(PS_row,var+1);
        thrust(resp_row,var)       = A(PS_row,num_vars+2);
        sfc(resp_row,var) = A(PS_row,num_vars+3);
        massFlow(resp_row,var) = A(PS_row,num_vars+4);
        thermalEff(resp_row,var) = A(PS_row,num_vars+5);
        PS_row = PS_row+1;
    end
end

% ----------------------
% ----------------------
% PLOTTING
% ----------------------
% ----------------------

steps(:,1) = -num_steps:num_steps;
input_vars = {'bypass','fanPstag','fanEff','compressEff','compressPratio','burnerPstag','burnerEff','turbineEffPoly','turbineEffShaft','Abypass2Acore'};
lower_bound = [0.57     2.91       0.82     0.84          24.0             0.92          0.94        0.83             0.95              0.15];
upper_bound = [0.63     3.21       0.86     0.9           25.0             0.98          0.99        0.89             0.99              0.4];

for ii = 1:size(thrust,1)
    for jj = 1:size(thrust,2)
        if(min(Vars(ii,jj) - lower_bound(jj)) > 0 && max(Vars(ii,jj) - upper_bound(jj)) < 0)
            thrust_inbounds(ii,jj) = thrust(ii,jj);
            sfc_inbounds(ii,jj) = sfc(ii,jj);
            massFlow_inbounds(ii,jj) = massFlow(ii,jj);
            thermalEff_inbounds(ii,jj) = thermalEff(ii,jj);
        else
            thrust_inbounds(ii,jj) = NaN;
            sfc_inbounds(ii,jj) = NaN;
            massFlow_inbounds(ii,jj) = NaN;
            thermalEff_inbounds(ii,jj) = NaN;
        end
    end
end

figure;
hold on;
% cstr = {'r^-', 'gv-', 'b>-', 'k.-', 'ms-', 'yd-', 'c--', 'r+-', 'g*-', 'mp-'};
cstr = {'r', 'g', 'b', 'k', 'm', 'y', 'c', 'r--', 'g--', 'm--'};
for ii = 1:num_vars
%     plot(steps, thrust(:,ii),cstr{ii});
    plot(steps, thrust_inbounds(:,ii),cstr{ii},'linewidth',3);
end
xlabel('Parameter step');
ylabel('Thurst');
legend(input_vars);
title('Centered parameter study for Thrust');
print('-dpdf','Thrust-centered-param-study-overlay.pdf');

figure;
hold on;
for ii = 1:num_vars
    plot(steps, sfc_inbounds(:,ii),cstr{ii},'linewidth',3);
end
xlabel('Parameter step');
ylabel('sfc');
legend(input_vars);
title('Centered parameter study for sfc');
print('-dpdf','sfc-centered-param-study-overlay.pdf');

figure;
hold on;
for ii = 1:num_vars
    plot(steps, massFlow_inbounds(:,ii),cstr{ii},'linewidth',3);
end
xlabel('Parameter step');
ylabel('massFlowRate');
legend(input_vars);
title('Centered parameter study for massFlowRate');
print('-dpdf','massFlowRate-centered-param-study-overlay.pdf');

figure;
hold on;
for ii = 1:num_vars
    plot(steps, thermalEff_inbounds(:,ii),cstr{ii},'linewidth',3);
end
xlabel('Parameter step');
ylabel('thermalEfficiency');
legend(input_vars);
title('Centered parameter study for thermalEfficiency');
print('-dpdf','thermalEfficiency-centered-param-study-overlay.pdf');
