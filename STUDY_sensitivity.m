%% 
% Initial subspace / sensitivity study for turbofan
clear all; close all;

%%
m = 10;
M = 100;
X = 2*rand(M,m)-1;

thrusts = zeros(M,1);
sfcs = zeros(M,1);
efficiencys = zeros(M,1);

for i=1:M
    
    [thrust, sfc, thermalEfficiency] = wrapperTurbofan(X(i,:));
    thrusts(i) = thrust.total;
    sfcs(i) = sfc;
    efficiencys(i) = thermalEfficiency;
    
end


%%
% thrust
A = [ones(M,1) X];
u = A\thrusts;
w = u(2:end)/norm(u(2:end));

figure(1);
plot(X*w,thrusts,'ko','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','r');
axis square; grid on;
set(gca,'FontSize',14);
xlabel('Active variable');
ylabel('Thrust');
print('figs/sp_thrust','-dpdf','-r300');

figure(2);
plot(1:m,w,'ko','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','b');
axis square; grid on;
set(gca,'FontSize',14);
xlabel('Parameter index');
ylabel('Weights (Thrust)');
ylim([-1 1]);
print('figs/weights_thrust','-dpdf','-r300');

%%
% sfc
A = [ones(M,1) X];
u = A\sfcs;
w_sfc = u(2:end)/norm(u(2:end));

figure(3);
plot(X*w_sfc,sfcs,'ko','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','r');
axis square; grid on;
set(gca,'FontSize',14);
xlabel('Active variable');
ylabel('sfc');
print('figs/sp_sfc','-dpdf','-r300');

figure(4);
plot(1:m,w_sfc,'ko','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','b');
axis square; grid on;
set(gca,'FontSize',14);
xlabel('Parameter index');
ylabel('Weights (sfc)');
ylim([-1 1]);
print('figs/weights_sfc','-dpdf','-r300');

%%
% efficiency
A = [ones(M,1) X];
u = A\efficiencys;
w_eff = u(2:end)/norm(u(2:end));

figure(5);
plot(X*w_eff,efficiencys,'ko','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','r');
axis square; grid on;
set(gca,'FontSize',14);
xlabel('Active variable');
ylabel('Thermal Efficiency');
print('figs/sp_efficiency','-dpdf','-r300');

figure(6);
plot(1:m,w_eff,'ko','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','b');
axis square; grid on;
set(gca,'FontSize',14);
xlabel('Parameter index');
ylabel('Weights (Thermal Efficiency)');
ylim([-1 1]);
print('figs/weights_efficiency','-dpdf','-r300');


