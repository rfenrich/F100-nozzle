clc; clear all; close all;


path = './';
file_list = {'Turbofan-crossvalid-L1.tab','Turbofan-crossvalid-L2.tab','Turbofan-crossvalid-L3.tab','Turbofan-crossvalid-L4.tab','Turbofan-crossvalid-L5.tab'};
path2 = './output/';
file_list2 = {'TurbofanUQ-crossvalid-L1.out','TurbofanUQ-crossvalid-L2.out','TurbofanUQ-crossvalid-L3.out','TurbofanUQ-crossvalid-L4.out','TurbofanUQ-crossvalid-L5.out'};

response_functions = {'thrust', 'sfc', 'massFlowRate', 'thermalEfficiency'};

A = importdata(sprintf('%s%s',path,'Turbofan-samples.tab.bak'));
data_truth = A.data;

Nd = 10;
Nrf = 4;

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
        err(nn,ii) = norm(diff_data(:,Nd+ii));
        rel_err(nn,ii) = err(nn,ii)/mean(data_truth(:,Nd+ii));
    end
end


%%
% Plots

figure;
semilogx(N_coeff,rel_err*100);
xlabel('# Expansion Terms');
ylabel('% Error');
legend(response_functions);
print('-dpdf','plots/Rel_err.pdf');


figure;
for nn = 1:numel(response_functions)
    subplot(2,2,nn);
    semilogx(N_coeff,err(:,nn));
    xlabel('# Expansion Terms');
    ylabel('$\vert\vert Error\vert\vert_2$','interpreter','latex');
    title(response_functions{nn});
end
print('-dpdf','plots/Abs_err.pdf');
