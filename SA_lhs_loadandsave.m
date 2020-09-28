clear all



%% Load data
filePath = dir(['*iter*', filesep, '*iter*.mat']);

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



%% Process data
% Get unique input parameters
params_fields = {'a', 'b', 'r0', 'mu'};
params_mat = [[params.a]', [params.b]', [params.r0]', [params.mu]'];
[params_unique, ~, params_index] = unique(params_mat, 'rows');
nPoint = size(params_unique, 1);

% For data fields of interest mean/std/sem output responses for each unique
% parameter point
data_fields = {'V_peakTime', 'V_peak', 't_end', 'statR', 'maxmaxR', 'statDiv', 'statD'};
nField = length(data_fields);
data_mat = zeros(nParams, nField);
data_mean = zeros(nPoint, nField);
data_std = zeros(nPoint, nField);
data_sem = zeros(nPoint, nField);
for iField=1:length(data_fields)
    data_mat(:,iField) = [data.(data_fields{iField})];
    data_mean(:,iField) = accumarray(params_index, data(:,iField), [], @mean);
    data_std(:,iField) = accumarray(params_index, data(:,iField), [], @std);
    data_sem(:,iField) = accumarray(params_index, data(:,iField), [], @(x)std(x)/sqrt(length(x)));
end
          
% PRCC
[PRCC_mat, pval_mat] = partialcorri(data_mean, params_unique, 'type', 'Spearman');



%% Save
save('SA_lhs.mat', ...
     'params_mat', 'params_index', 'params_fields', 'data_mat', 'data_fields', ... % data per simulation
     'params_unique', 'data_mean', 'data_std', 'data_sem', ... % data per point in hypercube
     'PRCC_mat', 'pval_mat') % global sensitivity analysis results