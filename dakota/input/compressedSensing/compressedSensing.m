clc; clear all; close all;

N_M = 10;
N_d = 10;
N_tot = 10000;
N_T = 2000;
M_all = round(logspace(log10(20),3,N_M));
N_k = 10;

% N_tot = 1000;
% N_T = 400;
% N_M = 10;
% N_d = 10;
% M_all = round(logspace(log10(10),log10(500),N_M));
% N_k = 1;

response_ind = [2];
response_func = {'thrust', 'sfc', 'massFlowRate', 'thermalEfficiency', 'turbineTstag', 'fanMach'};

%%
run_all = false;
% run_all = true;
if(run_all)
    % Seed random number generator for repeatable results:
    rng(1);
    
    X = rand(N_d,N_tot);
    
    X_sigma = [0.01, 0.05, 0.01, 0.01, 0.1, 0.01, 0.01, 0.015, 0.01, 0.007];
    X_mu = [0.62, 3.06, 0.84, 0.87, 24.5, 0.95, 0.95, 0.85, 0.97, 0.4];
    X_l = [0.59, 2.99, 0.82, 0.84, 24.0, 0.92, 0.94, 0.83, 0.95, 0.38];
    X_u = [0.64, 3.14, 0.86, 0.88, 25.0, 0.98, 0.99, 0.87, 0.99, 0.42];
    
    X = X.*repmat(X_sigma',1,N_tot) + repmat(X_mu',1,N_tot);
    for ii = 1:N_d
        for jj = 1:N_tot
            if(X(ii,jj) < X_l(ii))
                X(ii,jj) = X_l(ii);
            elseif(X(ii,jj) > X_u(ii))
                X(ii,jj) = X_lu(ii);
            end
        end
    end
    
    % Write X to file:
    Xfilename = 'Xfile.dat';
    
    % Write S to file:
    dlmwrite(Xfilename,X','delimiter',' ','precision','%.16e');
    
    % Update input file:
    str = fileread('TurbofanCompressedSensing_X_all.in.template');
    str = strrep(str, '<X-points-filename>', sprintf('''%s''',Xfilename));
    
    fid = fopen('TurbofanCompressedSensing_X_all.in', 'w');
    fwrite(fid, str, '*char');
    fclose(fid);
    
    % Run dakota for all X:
    if(exist('dakota.rst.CompressedSensing', 'file') == 2)
        system('dakota -read_restart=dakota.rst.CompressedSensing -i TurbofanCompressedSensing_X_all.in -o ../../output/TurbofanCompressedSensing_X_all.out > /dev/null');
    else
        system('dakota -i TurbofanCompressedSensing_X_all.in -o ../../output/TurbofanCompressedSensing_X_all.out > /dev/null');
        % Save dakota restart file:
        system('cp dakota.rst dakota.rst.CompressedSensing');
    end
    
    template_filename = 'TurbofanCompressedSensing.in.template';
    base_filename = 'TurbofanCompressedSensing';
    
    
    X_data = importdata('TurbofanCompressedSensing_X_all.tab');
    X_data = X_data(:,[1:N_d N_d + response_ind]);
    
    for M = M_all
        min_diff = NaN;
        p = 1;
        for pp = 1:6
            Np = nchoosek(N_d + pp,pp);
            if(abs(Np - 10*M) < min_diff || isnan(min_diff))
                min_diff = abs(Np - 10*M);
                p = pp;
            end
        end
        
        for k = 1:N_k
            I_all = randperm(N_tot);
            I = I_all(1:M);
            
            S = X_data(I,:);
            X_temp = X;
            X_temp(:,I) = [];
            X_data_temp = X_data;
            X_data_temp(I,:) = [];
            I2_all = randperm(N_tot-M);
            I2 = I2_all(1:N_T);
            T = X_temp(:,I2);
            T_data = X_data_temp(I2,:);
            
            Sfilename = sprintf('Sfile_M%d_k%d.dat',M,k);
            Tfilename = sprintf('Tfile_M%d_k%d.dat',M,k);
            T_datafilename = sprintf('T_datafile_M%d_k%d.dat',M,k);
            inputFilename = sprintf('%s_M%d_k%d.in',base_filename,M,k);
            dakotaOutFilename = sprintf('%s_M%d_k%d.out',base_filename,M,k);
            outputFilename = sprintf('output_M%d_k%d.dat',M,k);
            
            % Write S to file:
            dlmwrite(Sfilename,S,'delimiter',' ','precision','%.16e');
            
            % Write T to file:
            dlmwrite(Tfilename,T','delimiter',' ','precision','%.16e');
            
            % Write T_data to file:
            dlmwrite(T_datafilename,T_data,'delimiter',' ','precision','%.16e');
            
            % Update input file:
            str = fileread(template_filename);
            str = strrep(str, '<p>', sprintf('%d',p));
            %         str = strrep(str, '<M>', sprintf('%d',M));
            str = strrep(str, '<S-points-filename>', sprintf('''%s''',Sfilename));
            str = strrep(str, '<T-points-filename>', sprintf('''%s''',Tfilename));
            str = strrep(str, '<output-filename>', sprintf('''%s''',outputFilename));
            str = strrep(str, '<response>', sprintf('''%s''',response_func{response_ind}));
            
            fid = fopen(inputFilename, 'w');
            fwrite(fid, str, '*char');
            fclose(fid);
            
            % Run dakota:
            disp(sprintf('dakota -read_restart=dakota.rst.CompressedSensing -i %s -o ../../output/%s > /dev/null',inputFilename,dakotaOutFilename));
            tic;
            %         system(sprintf('dakota -read_restart=dakota.rst.CompressedSensing -i %s -o ../../output/%s > /dev/null',inputFilename,dakotaOutFilename));
            system(sprintf('dakota -i %s -o ../../output/%s > /dev/null',inputFilename,dakotaOutFilename));
            toc;
        end
    end
end

%%
% Post-process data:

nn = 0;
for M = M_all
    nn = nn + 1;
    min_diff = NaN;
    p = 1;
    for pp = 1:6
        Np_temp = nchoosek(N_d + pp,pp);
        if(abs(Np_temp - 10*M) < min_diff || isnan(min_diff))
            min_diff = abs(Np_temp - 10*M);
            Np = Np_temp;
            p = pp;
        end
    end
    
    for k = 1:N_k
        T_datafilename = sprintf('T_datafile_M%d_k%d.dat',M,k);
        outputFilename = sprintf('output_M%d_k%d.dat',M,k);
        
        if(exist(outputFilename, 'file') == 2)
            data_truth = importdata(T_datafilename);
        else
            err(nn,k,:) = NaN;
            break;
        end
        
        if(exist(outputFilename, 'file') == 2)
            A = importdata(outputFilename);
            if(isstruct(A))
                data_surrogate = A.data;
            else
                err(nn,k,:) = NaN;
                break;
            end
        else
            err(nn,k,:) = NaN;
            break;
        end
        
        
        diff_data = data_surrogate - data_truth;
        
        Nrf = size(diff_data,2) - N_d;
        Nsamp = size(diff_data,1);
        
        for ii = 1:Nrf
            err(nn,k,ii) = norm(diff_data(:,N_d+ii))/sqrt(Nsamp);
        end
    end
    
    N_coeff(nn) = Np;
end

mean_err = mean(err,2);
std_err = std(err,0,2);

med_err = zeros(size(err,1),size(err,3));
Q1_err = zeros(size(err,1),size(err,3));
Q2_err = zeros(size(err,1),size(err,3));
for ii = 1:size(err,1)
    for jj = 1:size(err,3)
        err_sort = sort(err,2);
        med_err(ii,jj) = squeeze(median(err_sort(ii,:,jj)));
        Q1_err(ii,jj) = median(squeeze(err_sort(ii,err_sort(ii,:,jj)<med_err(ii,jj),jj)));
        Q2_err(ii,jj) = median(squeeze(err_sort(ii,err_sort(ii,:,jj)>med_err(ii,jj),jj)));
    end
end


figure;
% h1 = loglog(M_all,squeeze(err(:,:,1)),'bx');
h1 = loglog(M_all,med_err(:,1),'b');
hold on;
loglog(M_all,Q1_err(:,1),'b--');
loglog(M_all,Q2_err(:,1),'b--');
% h = errorbar(N_coeff,mean_err(:,1), std_err(:,1));
% h = errorbar(M_all,mean_err(:,1), std_err(:,1));
% set(get(h,'Parent'), 'YScale', 'log');
% set(get(h,'Parent'), 'XScale', 'log');

xlabel('Model Evaluations');
ylabel('$\frac{\vert\vert Error\vert\vert_2}{\sqrt{N}}$','interpreter','latex');
title('Specific Fuel Consumption');




path = '../../output/';
file_list = {'Turbofan-crossvalid-L1.tab','Turbofan-crossvalid-L2.tab','Turbofan-crossvalid-L3.tab','Turbofan-crossvalid-L4.tab'};%,'Turbofan-crossvalid-L5.tab'};
path2 = '../../output/';
file_list2 = {'TurbofanUQ-crossvalid-L1.out','TurbofanUQ-crossvalid-L2.out','TurbofanUQ-crossvalid-L3.out','TurbofanUQ-crossvalid-L4.out'};%,'TurbofanUQ-crossvalid-L5.out'};

response_functions = {'thrust', 'sfc', 'massFlowRate', 'thermalEfficiency'};

A = importdata(sprintf('%s%s',path,'Turbofan-samples.tab.bak'));
data_truth = A.data;

Nd = 10;
Nrf = 4;
Nsamp = 2000;

N_coeff2 = zeros(1,numel(file_list2));
for nn = 1:numel(file_list2)
    [stat out] = system(sprintf('grep ''Total number'' %s%s',path2,file_list2{nn}));
    temp = textscan(out,'%s');
    N_coeff2(nn) = str2num(char(temp{1}(6)));
end

for nn = 1:numel(file_list)
    A = importdata(sprintf('%s%s',path,file_list{nn}));
    data_surrogate = A.data;
    
    
    diff_data = data_surrogate - data_truth;
    
    for ii = 1:Nrf
        err2(nn,ii) = norm(diff_data(:,Nd+ii))/sqrt(Nsamp);
    end
end

hold on;
nn = 2;
h2 = loglog(N_coeff2,err2(1:4,nn),'r');
legend([h1,h2],{'Adaptive Basis, Compressed Sensing','Sparse Grid'});
matlab2tikz('compressedSensingComparison.tex','standalone',true);