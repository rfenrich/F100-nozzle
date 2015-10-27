clc; clear all; close all;

response_functions = {'thrust', 'sfc', 'massFlowRate', 'thermalEfficiency'};
lit_values = [64900, 0.73, 102, NaN];
% file_list = {'TurbofanUQ-pce-L0.out','TurbofanUQ-pce-L1.out','TurbofanUQ-pce-L2.out','TurbofanUQ-pce-L3.out','TurbofanUQ-pce-L4.out','TurbofanUQ-pce-L5.out'};
file_list = {'TurbofanUQ-pce-L0.out','TurbofanUQ-pce-L1.out','TurbofanUQ-pce-L2.out','TurbofanUQ-pce-L3.out','TurbofanUQ-pce-L4.out'};
levels = 0:1:numel(file_list)-1;
input_vars = {'bypass','fanPstag','fanEff','compressEff','compressPratio','burnerPstag','burnerEff','turbineEffPoly','turbineEffShaft','Abypass2Acore'};

path = './output/';


response_mean = zeros(numel(response_functions),numel(file_list));
response_std = zeros(numel(response_functions),numel(file_list));
sobol_main_ind = zeros(numel(response_functions),numel(file_list),numel(input_vars));
sobol_total_ind = zeros(numel(response_functions),numel(file_list),numel(input_vars));


N_coeff = zeros(1,numel(file_list));
for nn = 1:numel(file_list)
    [stat out] = system(sprintf('grep ''Total number'' %s%s',path,file_list{nn}));
    temp = textscan(out,'%s');
    N_coeff(nn) = str2num(char(temp{1}(6)));
end


alpha_len = zeros(numel(file_list),numel(response_functions));
alpha_raw = 0;
index_raw = zeros(numel(file_list), numel(response_functions), max(N_coeff), numel(input_vars));
max_ind_coeff = numel(file_list);
for nn = 1:numel(file_list)
    fileID = fopen(sprintf('%s%s',path,file_list{nn}),'r');
    
    if(nn <= max_ind_coeff)
        for ii = 1:numel(response_functions)
            % Search for normalized coefficients:
            while (~feof(fileID))
                line = fgetl(fileID);
                if(~isempty(findstr(line,'Normalized coefficients of Polynomial Chaos Expansion for')))
                    % found it!
                    break;
                end
            end
            
            for jj = 1:numel(response_functions)
                if(~isempty(findstr(line,response_functions{jj})))
                    % Discard header:
                    line = fgetl(fileID);
                    line = fgetl(fileID);
                    
                    line = fgetl(fileID);
                    count = 1;
                    while(isempty(findstr(line,'Normalized coefficients')) && isempty(findstr(line,'--------------------')))
                        temp = textscan(line, '%s');
                        alpha(nn,jj,count) = str2num(char(temp{1}(1)));
                        
                        for mm = 2:numel(temp{1})
                            temp_str = char(temp{1}(mm));
                            index_raw(nn,jj,count,mm-1) = str2num(temp_str(2:end));
                        end
                        
                        line = fgetl(fileID);
                        count = count+1;
                    end
                    alpha_len(nn,jj) = count-1;
                end
            end
        end
    end
    fclose(fileID);
    
    
    fileID = fopen(sprintf('%s%s',path,file_list{nn}),'r');
    % Search for response function statistics:
    while (~feof(fileID))
        line = fgetl(fileID);
        if(~isempty(findstr(line,'Moment statistics for each response function:')))
            % found it!
            break;
        end
    end
    
    % discard header:
    line = fgetl(fileID);
    
    for ii = 1:numel(response_functions)
        line = fgetl(fileID);
        for jj = 1:numel(response_functions)
            if(~isempty(findstr(line,response_functions{jj})))
                line = fgetl(fileID);
                mean_std_tmp = textscan(line, ' expansion: %f %f');
                response_mean(jj,nn) = mean_std_tmp{1};
                response_std(jj,nn) = mean_std_tmp{2};
                
                line = fgetl(fileID);
            end
        end
    end
    
    
    % Search for global sensitivities
    while (~feof(fileID))
        line = fgetl(fileID);
        if(~isempty(findstr(line,'Global sensitivity indices for each response function:')))
            % found it!
            break;
        end
    end
    
    if(norm(response_std(:,nn)) > eps)
        for ii = 1:numel(response_functions)
            line = fgetl(fileID);
            for jj = 1:numel(response_functions)
                if(~isempty(findstr(line,response_functions{jj})))
                    % discard header:
                    line = fgetl(fileID);
                    
                    for mm = 1:numel(input_vars)
                        line = fgetl(fileID);
                        temp = textscan(line, '%f %f %s');
                        for kk = 1:numel(input_vars)
                            if(~isempty(findstr(char(temp{3}),input_vars{kk})))
                                sobol_main_ind(jj,nn,kk) = temp{1};
                                sobol_total_ind(jj,nn,kk) = temp{2};
                            end
                        end
                    end
                    
                    while (~feof(fileID) && jj+1 <= numel(response_functions))
                        line = fgetl(fileID);
                        if(~isempty(findstr(line,response_functions{jj+1})))
                            % found it!
                            break;
                        end
                    end
                end
            end
        end
    end
    
    
    fclose(fileID);
end

%%
% Plots

figure;
for ind = 1:numel(response_functions)
    subplot(2,2,ind);
    errorbar(levels,response_mean(ind,:), response_std(ind,:));
    hold on;
    plot(levels,lit_values(ind)*ones(1,numel(file_list)),'r');
    xlabel('Sparse Grid Levels');
    ylabel(sprintf('%s',response_functions{ind}));
end

for ind1 = 2:numel(file_list)
    figure;
    for ind2 = 1:numel(response_functions)
        subplot(2,2,ind2);
        bar(squeeze(sobol_main_ind(ind2,ind1,:)));
        title(response_functions{ind2});
%         set(gca, 'XTickLabel',input_vars, 'XTick',1:numel(input_vars));
        xticklabel_rotate([1:numel(input_vars)],90,input_vars,'interpreter','none');
        ylabel('Main Sobol'' Indices');
    end
    if(ind1 == 5)
        print('-dpdf',  'plots/MainSobolInd.pdf');
    end
end

for ind1 = 2:numel(file_list)
    figure;
    for ind2 = 1:numel(response_functions)
        subplot(2,2,ind2);
        bar(squeeze(sobol_total_ind(ind2,ind1,:)));
        title(response_functions{ind2});
%         set(gca, 'XTickLabel',input_vars, 'XTick',1:numel(input_vars));
        xticklabel_rotate([1:numel(input_vars)],90,input_vars,'interpreter','none');
        ylabel('Total Sobol'' Indices');
    end
    if(ind1 == 5)
%         print('-dpdf',  'plots/TotalSobolInd.pdf');
    end
    matlab2tikz(sprintf('plots/TotalSobolIndL-%d.tex', ind1-1),'standalone',true);
end
%%
error_mean = (response_mean(:,1:end-1)-repmat(response_mean(:,end), 1, size(response_mean,2)-1))./repmat(response_mean(:,end), 1, size(response_mean,2)-1);
figure;
semilogy(levels(1:end-1),abs(error_mean));
xlabel('Sparse Grid Levels');
ylabel('$\vert \% Error\vert$','interpreter','latex');
legend(response_functions);
title('Error of Response Function Mean');

error_std = (response_std(:,1:end-1)-repmat(response_std(:,end), 1, size(response_std,2)-1))./repmat(response_std(:,end), 1, size(response_std,2)-1);
figure;
semilogy(levels(1:end-1),abs(error_std));
xlabel('Sparse Grid Levels');
ylabel('$\vert \% Error\vert$','interpreter','latex');
legend(response_functions);
title('Error of Response Function Std. Dev.');

error_sobol_main = sobol_main_ind(:,1:end-1,:) - repmat(sobol_main_ind(:,end,:),[ 1, size(sobol_main_ind,2)-1,1]);
figure;
for ind1 = 1:numel(response_functions)
    subplot(2,2,ind1);
    semilogy(levels(1:end-1),abs(squeeze(error_sobol_main(ind1,:,:))));
    title(response_functions{ind1});
    xlabel('Sparse Grid Levels');
    ylabel('$\vert Error\vert$ -- Main Sobol'' Indices','interpreter','latex');
    if(ind1 == 1)
        legend(input_vars);
    end
end

error_sobol_total = sobol_total_ind(:,1:end-1,:) - repmat(sobol_total_ind(:,end,:),[ 1, size(sobol_total_ind,2)-1,1]);
figure;
for ind1 = 1:numel(response_functions)
    subplot(2,2,ind1);
    semilogy(levels(1:end-1),abs(squeeze(error_sobol_total(ind1,:,:))));
    title(response_functions{ind1});
    xlabel('Sparse Grid Levels');
    ylabel('$\vert Error\vert$ -- Total Sobol'' Indices','interpreter','latex');
    if(ind1 == 1)
        legend(input_vars);
    end
end


%%
% Plot normalized PCE coefficients
total_expansion = importdata('ref_eo16_d10.mat');

for ii = 1:16
    total_ind_levels(ii) = factorial(numel(input_vars) + ii)/(factorial(numel(input_vars))*factorial(ii));
end

tic;
for nFile = 2:max_ind_coeff
    figure;
    for nFunc = 1:numel(response_functions)
        index_raw_sub = squeeze(index_raw(nFile,nFunc,:,:));
        alpha_sub = squeeze(alpha(nFile,nFunc,1:alpha_len(nFile,nFunc)));
        index_alpha = zeros(alpha_len(nFile,nFunc),1);
        for ii=1:alpha_len(nFile,nFunc)
            for jj=1:size(total_expansion, 1)
                if(isequal(index_raw_sub(ii,:), total_expansion(jj,:)))
                    index_alpha(ii) = jj;
                    break;
                end
            end
            if(index_alpha(ii) == 0)
                index_alpha(ii) = NaN;
            end
        end
        
        subplot(2,2,nFunc);
        loglog(index_alpha, abs(alpha_sub),'x');
        hold on;
        for ii = 1:16
            loglog([total_ind_levels(ii) total_ind_levels(ii)], [min(abs(alpha_sub)) max(abs(alpha_sub))],'k--');
        end
        title(sprintf('Sparse Grid L-%d, %s', nFile-1, response_functions{nFunc}));
        xlabel('Basis id (total-order basis)');
        ylabel('$\vert \alpha_i\vert$','interpreter','latex');
    end
%     print('-dpdf',  sprintf('plots/SparseGridL-%d.pdf', nFile-1));
    matlab2tikz(sprintf('plots/SparseGridL-%d.tex', nFile-1),'standalone',true);
end
toc;

% for nFile = 2:max_ind_coeff
%     figure;
%     for nFunc = 1:numel(response_functions)
%         index_raw_sub = squeeze(index_raw(nFile,nFunc,:,:));
%         alpha_sub = squeeze(alpha(nFile,nFunc,1:alpha_len(nFile,nFunc)));
%         index_alpha_simple = zeros(alpha_len(nFile,nFunc),1);
%         for ii=1:alpha_len(nFile,nFunc)
%             index_alpha_simple(ii) = sum(index_raw_sub(ii,:));
%         end
%         
%         subplot(2,2,nFunc);
%         semilogy(index_alpha_simple, abs(alpha_sub),'x');
%         title(sprintf('Sparse Grid L-%d, %s', nFile-1, response_functions{nFunc}));
%         xlabel('Basis id (total-order basis)');
%         ylabel('$\vert \alpha_i\vert$','interpreter','latex');
%     end
%     print('-dpdf',  sprintf('SparseGridL-%d-Simple.pdf', nFile-1));
% end


