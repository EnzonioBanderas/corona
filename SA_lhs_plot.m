clear all

% natural variance depending on position in hypercube?
% increase number of points in hypercube to improve infection rate
% significance?
% 
%% Define
tempfig = figure('units','normalized','outerposition',[0 0 1 1]);
cmap = colormap(lines);
close(tempfig)

set(groot,'defaulttextinterpreter','none');  
set(groot,'defaultAxesTickLabelInterpreter','none');  
set(groot,'defaultLegendInterpreter','none');

% loop over parameter varies here?
parameter_varies = {'lhs_lin2_U0_1e4_nPoint_1e3_nK_'};
% for iPV = 1:length(parameter_varies)
% parameter_vary = 'a_lin'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iPV = 1;
parameter_vary = parameter_varies{iPV};
parameter_vary_split = strsplit(parameter_vary, '_');



%% Load data
filePath = dir(['D:\Users\enzo\Downloads\*', parameter_vary, '*\*\*iter3_*.mat']);
nFile = length(filePath);
data_cell = cell(1, nFile);
params_cell = cell(1, nFile);
ID_cell = cell(1, nFile);
nIter = zeros(1, nFile);
for iFile = 1:nFile
  load(fullfile(filePath(iFile).folder, filePath(iFile).name))
  data_cell{iFile} = data;
  nIter(iFile) = length(data);
  ID_cell{iFile} = repmat(iFile, [1, nIter(iFile)]);
  params_cell{iFile} = repmat(params, [1, nIter(iFile)]);
  for iIter = 1:nIter(iFile)
      data_cell{iFile}(iIter).maxmaxR = max(data_cell{iFile}(iIter).maxR);
  end
  clear data params
end
data = [data_cell{:}];
params = [params_cell{:}];
ID = [ID_cell{:}];
nParams = length(params);

% stationary replication rate
% maximum replication rate
% stationary evenness
% stationary hamming distance
% V peak time
% V peak value
% end of infection time

% what other measures are interesting?


% stationary fitness maximum
% stationary diversity hamming

%% Process data

% Get input parameter and sort
params_fields = {'a', 'b', 'r0', 'mu'};
params_mat = [[params.a]', [params.b]', [params.r0]', [params.mu]'];
input_mat = [[params.a]', [params.b]', [params.r0]', [params.mu]'];
[input, ~, accum_index] = unique(params_mat, 'rows');
nInput = size(input, 1);

% data_fields = fields(data);
data_fields = {'V_peakTime', 'V_peak', 't_end', 'statR', 'maxmaxR', 'statDiv', 'statD'};
nField = length(data_fields);
data_mean = zeros(nInput, nField);
data_std = zeros(nInput, nField);
data_sem = zeros(nInput, nField);
output_mat = zeros(nParams, nField);
for iField=1:length(data_fields)
    output_mat(:,iField) = [data.(data_fields{iField})];
    data_mean(:,iField) = accumarray(accum_index, [data.(data_fields{iField})], [], @mean);
    data_std(:,iField) = accumarray(accum_index, [data.(data_fields{iField})], [], @std);
    data_sem(:,iField) = accumarray(accum_index, [data.(data_fields{iField})], [], @(x)std(x)/sqrt(length(x)));
end

% % Filter out nan values which result due to certain bins missing values (bins can not contain both nan values and valid values right?)
% input_bin_logical = ~isnan(input_bin);
% input_bin = input_bin(input_bin_logical,:);
% data_mean = data_mean(input_bin_logical,:);
% data_std = data_std(input_bin_logical,:);
% data_sem = data_sem(input_bin_logical,:);
          
% PRCC
[PRCC_mat, pval_mat] = partialcorri(data_mean, input);


% [PRCC_mat, pval_mat] = partialcorri(output_mat, input_mat);



%% Plot PRCC heatmap
fig1 = figure();
movegui(fig1,'west');

PRCC_mat_plot = PRCC_mat;
PRCC_mat_plot(~(pval_mat<0.05)) = nan;
hm = heatmap(PRCC_mat_plot);
hm.XDisplayLabels = params_fields;
hm.XLabel = 'Input parameters';
hm.YDisplayLabels = {'V peakTime', 'V peak', 't end', 'statR', 'maxmaxR', 'statDiv', 'statD'};
hm.YLabel = 'Output responses';
hm.Title = 'PRCC matrix';
colormap(parula)



%% Plot pval heatmap
fig2 = figure();
movegui(fig2,'east');

hm = heatmap(pval_mat);
hm.XDisplayLabels = params_fields;
hm.XLabel = 'Input parameters';
hm.YDisplayLabels = {'V peakTime', 'V peak', 't end', 'statR', 'maxmaxR', 'statDiv', 'statD'};
hm.YLabel = 'Output responses';
hm.Title = 'pval matrix';
colormap(parula)



%% Plot barplots
% fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
figure('units','normalized','outerposition',[0 0 1 1])
sgtitle('PRCC values')

for iP = 1:length(params_fields)
    subplot(2,2,iP)
    bp(iP) = bar(PRCC_mat(:,iP));
    xticklabels(data_fields)
%     xticks(1:1:length(data_fields))
    xtickangle(45)
    YLIM = max(abs(PRCC_mat(:,iP)));
%     ylim([-YLIM, YLIM])
    ylim([-1, 1])
    title(params_fields{iP})
    ylabel('PRCC')
    xlabel('Output responses')
end




% %% Plot input output distribution
% fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
% 
% N = 50;
% subplot(2,1,1)
% % histogram(log10(input), N)
% % xlabel('Input parameter')
% histogram(input, N, 'Normalization', 'pdf')
% title([parameter_vary,' input'], 'Interpreter', 'none')
% xlabel('Input parameter')
% % xlabel('Logarithmic input')
% % xlabel('Logarithmic infection rate (log(a))')
% ylabel('Probability density function')
% % set(gca, 'XScale', 'log')
% subplot(2,1,2)
% hold on
% for iField=1:length(data_fields)
%     data_perfield = [data.(data_fields{iField})];
%     histogram(data_perfield/max(data_perfield), N, 'Normalization', 'pdf')
% end
% % histogram([data.V_peakTime], N), hold on
% % histogram([data.t_end], N)
% title([parameter_vary,' output'], 'Interpreter', 'none')
% xlabel('Normalized output response')
% ylabel('Probability density function')
% % output_str = {'V peak time', ...
% %               't end'};
% legend(legend_str, 'Interpreter', 'none', 'Location', 'bestoutside')
% 
% saveas(fig1, ['Output/SA/',parameter_vary,'_inputOutput'])
% 
% subplot(2,1,1)
% legend({'Input'}, 'Interpreter', 'none', 'Location', 'bestoutside')
% subplot(2,1,2)
% legend(legend_str, 'Interpreter', 'none', 'Location', 'bestoutside')
% saveas(fig1, ['Output/SA/',parameter_vary,'_inputOutput.png'])
% 
% 
% %% Plot output ratio to max (y) over input (x)
% fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
% 
% % out = cell(1,length(output_cell)); 
% for i=1:length(legend_str)
%     maxratioLines(i) = plot_line_shaded(input, ...
%         data_mean(:,i)/max(data_mean(:,i), [], 1), ...
%         data_sem(:,i)/max(data_mean(:,i), [], 1), ...
%         cmap(i,:));
%     hold on
% %     plot(input, data_mean_sorted./max(data_mean_sorted, [], 1)), hold on
% end
% 
% if strcmp(parameter_vary_split{2}, 'lin')
% else
%     set(gca, 'XScale', 'log')
% end
% legend(maxratioLines, legend_str, 'Interpreter', 'none')
% xlabel('Input parameter')
% ylabel('Normalized output response')
% title(parameter_vary, 'Interpreter', 'none')
% xlim([min(input), max(input)])
% 
% saveas(fig2, ['Output/SA/',parameter_vary,'_normalizedResponse'])
% 
% legend(maxratioLines, legend_str, 'Interpreter', 'none', 'Location', 'bestoutside')
% saveas(fig2, ['Output/SA/',parameter_vary,'_normalizedResponse.png'])
% 
% 
% %% Plot V peak time and infection time over input
% fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
% 
% count = 0;
% for i=[1,3]
%     count = count+1;
%     plotted_lines(count) = plot_line_shaded(input, ...
%     data_mean(:,i), data_sem(:,i), cmap(i,:)); 
%     hold on
% end
% 
% legend(plotted_lines, {'V peak time', 't end'})
% xlabel('Input parameter')
% ylabel('Time (days)')
% title(parameter_vary, 'Interpreter', 'none')
% 
% if strcmp(parameter_vary_split{2}, 'lin')
% else
%     set(gca, 'XScale', 'log')
% end
% xlim([min(input), max(input)])
% 
% saveas(fig3, ['Output/SA/',parameter_vary,'_timeResponse'])
% 
% legend(plotted_lines, {'V peak time', 't end'}, 'Interpreter', 'none', 'Location', 'bestoutside')
% saveas(fig3, ['Output/SA/',parameter_vary,'_timeResponse.png'])
% 
% 
% % end
