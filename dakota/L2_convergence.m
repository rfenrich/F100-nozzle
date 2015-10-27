% clc; clear all; close all;


path = './output/';
file_list = {'Turbofan-crossvalid-L1.tab','Turbofan-crossvalid-L2.tab','Turbofan-crossvalid-L3.tab','Turbofan-crossvalid-L4.tab'};%,'Turbofan-crossvalid-L5.tab'};
path2 = './output/';
file_list2 = {'TurbofanUQ-crossvalid-L1.out','TurbofanUQ-crossvalid-L2.out','TurbofanUQ-crossvalid-L3.out','TurbofanUQ-crossvalid-L4.out'};%,'TurbofanUQ-crossvalid-L5.out'};

response_functions = {'thrust', 'sfc', 'massFlowRate', 'thermalEfficiency'};

A = importdata(sprintf('%s%s',path,'Turbofan-samples.tab.bak'));
data_truth = A.data;

Nd = 10;
Nrf = 4;
Nsamp = 2000;

N_coeff = zeros(1,numel(file_list2));
for nn = 1:numel(file_list2)
    [stat out] = system(sprintf('grep ''Total number'' %s%s',path2,file_list2{nn}));
    temp = textscan(out,'%s');
    N_coeff(nn) = str2num(char(temp{1}(6)));
end

for nn = 1:numel(file_list)
    A = importdata(sprintf('%s%s',path,file_list{nn}));
    data_surrogate = A.data;
    
    
    diff_data = data_surrogate - data_truth;
    
    for ii = 1:Nrf
        err(nn,ii) = norm(diff_data(:,Nd+ii))/sqrt(Nsamp);
    end
end


%%
% Plots


% figure;
% for nn = 1:numel(response_functions)
%     subplot(2,2,nn);
%     loglog(N_coeff,err(1:4,nn));
%     xlabel('# Expansion Terms');
%     ylabel('$\sqrt{\frac{\vert\vert Error\vert\vert_2}{N}}$','interpreter','latex');
%     title(response_functions{nn});
% end
% % matlab2tikz('plots/L2-Convergence-Kickoff.tex','standalone',true);
% % print('-dpdf','plots/Abs_err.pdf');

%%
figure;
nn = 2;
loglog(N_coeff,err(1:4,nn));
xlabel('# Expansion Terms');
ylabel('$\sqrt{\frac{\vert\vert Error\vert\vert_2}{N}}$','interpreter','latex');
title(response_functions{nn});