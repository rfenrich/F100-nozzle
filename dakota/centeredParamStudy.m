clc; clear all;% close all;

response_functions = {'thrust', 'sfc', 'massFlowRate', 'thermalEfficiency'};
lit_values = [64900, 0.73, 102, NaN];
file_list = {'TurbofanUQ-centered-param-study.out'};
levels = 0:1:numel(file_list)-1;
input_vars = {'bypass','fanPstag','fanEff','compressEff','compressPratio','burnerPstag','burnerEff','turbineEffPoly','turbineEffShaft','Abypass2Acore'};

lower_bound = [0.57     2.91       0.82     0.84          24.0             0.92          0.94        0.83             0.95              0.15];
upper_bound = [0.63     3.21       0.86     0.9           25.0             0.98          0.99        0.89             0.99              0.4];

for nn = 1:numel(file_list)
    fileID = fopen(file_list{nn},'r');
    
    % Search for evaluation parameters:
    while (~feof(fileID))
        line = fgetl(fileID);
        if(~isempty(findstr(line,'Begin Evaluation')))
            temp = textscan(line, '%s');
            eval_num = str2num(char(temp{1}(3)));
            
            % discard header:
            line = fgetl(fileID);
            line = fgetl(fileID);
            
            for ii = 1:numel(input_vars)
                line = fgetl(fileID);
                for jj = 1:numel(input_vars)
                    if(~isempty(findstr(line,input_vars{jj})))
                        if(jj == 1 && ~isempty(findstr(line,input_vars{end})))
                            continue;
                        end
                        temp = textscan(line, '%s');
                        param(eval_num,jj) = str2num(char(temp{1}(1)));
                        param_inbounds(eval_num,jj) = NaN;
                    end
                end
            end
            
            % Check bounds for valid paramaters:
            if(min(param(eval_num,:) - lower_bound) > 0 && max(param(eval_num,:) - upper_bound) < 0)
                % Store
                param_inbounds(eval_num,:) = param(eval_num,:);
            end
        end
        
        if(~isempty(findstr(line,'Active response data for evaluation')))
            temp = textscan(line, '%s');
            temp = char(temp{1}(6));
            temp = temp(1:end-1);
            eval_num = str2num(temp);
            
            % discard header:
            line = fgetl(fileID);
            
            for ii = 1:numel(response_functions)
                line = fgetl(fileID);
                for jj = 1:numel(response_functions)
                    if(~isempty(findstr(line,response_functions{jj})))
                        temp = textscan(line, '%s');
                        response_val(eval_num,jj) = str2num(char(temp{1}(1)));
                        response_val_inbounds(eval_num,jj) = NaN;
                    end
                end
            end
            
            % Check bounds for valid paramaters:
            if(min(param(eval_num,:) - lower_bound) > 0 && max(param(eval_num,:) - upper_bound) < 0)
                % Store
                response_val_inbounds(eval_num,:) = response_val(eval_num,:);
            end
        end
    end
    
    
    fclose(fileID);
end

%%
% Plots

figure;
plot((param(:,:)-repmat(mean(param),size(param,1),1))./repmat(max((param(:,:)-repmat(mean(param),size(param,1),1))),size(param,1),1));

%%

for nn = 1:numel(response_functions)
    figure;
    plot(response_val(:,nn));
    hold on;
    plot(response_val_inbounds(:,nn),'k','linewidth',3);
    for ii = 1:size(response_val,1)/2000
        plot([ii*2000 + 1, ii*2000 + 1],[min(response_val(:,nn)) max(response_val(:,nn))],'k--');
        h = text(ii*2000 + 1 - 1500,min(response_val(:,nn)) + 0.05*(min(response_val(:,nn)) - max(response_val(:,nn))),input_vars{ii});
        set(h, 'rotation', 60);
    end
    title(response_functions{nn});
    xlabel('Evaluation Number');
    ylabel('Response Function');
    
    print('-dpdf',sprintf('%s-centered-param-study.pdf',response_functions{nn}));
end

