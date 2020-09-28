clear all



% %% Generate small test dataset
% N = 100;
% for x = 1:N
%     x_str = num2str(x);
%     Gillespie_func(x_str);
% end



%% Define
tempfig = figure('units','normalized','outerposition',[0 0 1 1]);
cmap = colormap(lines);
close(tempfig)

% loop over parameter varies here?
parameter_varies = {'a_lin', 'a_log', ...
    'b_lin', 'b_log', ...
    'c_lin', 'c_log', ...
    'r0_lin', 'r0_log', ...
    'mu_lin', 'mu_log', ...
    'V0_lin', 'V0_log'};
for iPV = 1:length(parameter_varies)
% parameter_vary = 'a_lin'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter_vary = parameter_varies{iPV};
parameter_vary_split = strsplit(parameter_vary, '_');
% if strcmp(parameter_vary_split{1}, 'c')
%     if strcmp(parameter_vary_split{2}, 'lin')
%         parameter_vary_split{2} = 'log';
%     else
%         parameter_vary_split{2} = 'lin';
%     end
% end


%% Load data
filePath = dir(['D:\Users\enzo\Downloads\*', parameter_vary, '\*\*_*.mat']);
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
      params_cell{iFile}(iIter).mu = params_cell{iFile}(iIter).mu(1);
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
switch parameter_vary_split{1}
    case 'a'
        input = [params.a];
    case 'b'
        input = [params.b];
    case 'c'
        input = [params.c];
    case 'r0'
        input = [params.r0];
    case 'mu'
        input = [params.mu];
    case 'V0'
        input = [params.V0];
end

% Bin data after sorting
nBin = round(length(ID)/(3*9));  % around 3*9 iterations per bin instead of 3
if strcmp(parameter_vary_split{2}, 'lin')
    [~, ~, ID_binned] = histcounts(input, nBin);
else
    [~, ~, ID_binned] = histcounts(log(input), nBin);
end
input_bin = accumarray(ID_binned', input, [], @median, nan);
% [input_bin_sorted, output_sort] = sort(input_bin);

% data_fields = fields(data);
data_fields = {'V_peakTime', 'V_peak', 't_end', 'statR', 'maxmaxR', 'statDiv', 'statD'};
nField = length(data_fields);
data_mean = zeros(nBin, nField); 
data_std = zeros(nBin, nField);
data_sem = zeros(nBin, nField);
for iField=1:length(data_fields)
    data_mean(:,iField) = accumarray(ID_binned', [data.(data_fields{iField})], [], @mean, nan);
    data_std(:,iField) = accumarray(ID_binned', [data.(data_fields{iField})], [], @std, nan);
    data_sem(:,iField) = accumarray(ID_binned', [data.(data_fields{iField})], [], @(x)std(x)/sqrt(length(x)), nan);
end

% Filter out nan values which result due to certain bins missing values (bins can not contain both nan values and valid values right?)
input_bin_logical = ~isnan(input_bin);
input_bin = input_bin(input_bin_logical,:);
data_mean = data_mean(input_bin_logical,:);
data_std = data_std(input_bin_logical,:);
data_sem = data_sem(input_bin_logical,:);

% Legend for normalized plots
legend_str = {['V_peakTime (', num2str(max(data_mean(:,1))), ')'], ...
              ['V_peak (', num2str(max(data_mean(:,2))), ')'], ...
              ['t_end (', num2str(max(data_mean(:,3))), ')'], ...
              ['statR (', num2str(max(data_mean(:,4))), ')'], ...
              ['maxmaxR (', num2str(max(data_mean(:,5))), ')'], ...
              ['statDiv (', num2str(max(data_mean(:,6))), ')'], ...
              ['statD (', num2str(max(data_mean(:,7))), ')']};
          
          

%% Plot input output distribution
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);

N = 50;
subplot(2,1,1)
% histogram(log10(input), N)
% xlabel('Input parameter')
histogram(input, N, 'Normalization', 'pdf')
title([parameter_vary,' input'], 'Interpreter', 'none')
xlabel('Input parameter')
% xlabel('Logarithmic input')
% xlabel('Logarithmic infection rate (log(a))')
ylabel('Probability density function')
% set(gca, 'XScale', 'log')
subplot(2,1,2)
hold on
for iField=1:length(data_fields)
    data_perfield = [data.(data_fields{iField})];
    histogram(data_perfield/max(data_perfield), N, 'Normalization', 'pdf')
end
% histogram([data.V_peakTime], N), hold on
% histogram([data.t_end], N)
title([parameter_vary,' output'], 'Interpreter', 'none')
xlabel('Normalized output response')
ylabel('Probability density function')
% output_str = {'V peak time', ...
%               't end'};
legend(legend_str, 'Interpreter', 'none', 'Location', 'bestoutside')

saveas(fig1, ['Output/SA/',parameter_vary,'_inputOutput'])

subplot(2,1,1)
legend({'Input'}, 'Interpreter', 'none', 'Location', 'bestoutside')
subplot(2,1,2)
legend(legend_str, 'Interpreter', 'none', 'Location', 'bestoutside')
saveas(fig1, ['Output/SA/',parameter_vary,'_inputOutput.png'])


%% Plot output ratio to max (y) over input (x)
fig2 = figure('units','normalized','outerposition',[0 0 1 1]);

% out = cell(1,length(output_cell)); 
for i=1:length(legend_str)
    maxratioLines(i) = plot_line_shaded(input_bin, ...
        data_mean(:,i)/max(data_mean(:,i), [], 1), ...
        data_sem(:,i)/max(data_mean(:,i), [], 1), ...
        cmap(i,:));
    hold on
%     plot(input, data_mean_sorted./max(data_mean_sorted, [], 1)), hold on
end

if strcmp(parameter_vary_split{2}, 'lin')
else
    set(gca, 'XScale', 'log')
end
legend(maxratioLines, legend_str, 'Interpreter', 'none')
xlabel('Input parameter')
ylabel('Normalized output response')
title(parameter_vary, 'Interpreter', 'none')
xlim([min(input_bin), max(input_bin)])

saveas(fig2, ['Output/SA/',parameter_vary,'_normalizedResponse'])

legend(maxratioLines, legend_str, 'Interpreter', 'none', 'Location', 'bestoutside')
saveas(fig2, ['Output/SA/',parameter_vary,'_normalizedResponse.png'])


%% Plot V peak time and infection time over input
fig3 = figure('units','normalized','outerposition',[0 0 1 1]);

count = 0;
for i=[1,3]
    count = count+1;
    plotted_lines(count) = plot_line_shaded(input_bin, ...
    data_mean(:,i), data_sem(:,i), cmap(i,:)); 
    hold on
end

legend(plotted_lines, {'V peak time', 't end'})
xlabel('Input parameter')
ylabel('Time (days)')
title(parameter_vary, 'Interpreter', 'none')

if strcmp(parameter_vary_split{2}, 'lin')
else
    set(gca, 'XScale', 'log')
end
xlim([min(input_bin), max(input_bin)])

saveas(fig3, ['Output/SA/',parameter_vary,'_timeResponse'])

legend(plotted_lines, {'V peak time', 't end'}, 'Interpreter', 'none', 'Location', 'bestoutside')
saveas(fig3, ['Output/SA/',parameter_vary,'_timeResponse.png'])


end
