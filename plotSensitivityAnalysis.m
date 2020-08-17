clear all

N = 100;
for x = 1:N
    x_str = num2str(x);
    Gillespie_func(x_str);
end

filePath = dir('data_iteration/data_iteration*');
nFile = length(filePath);
data_cell = cell(1, nFile);
params_cell = cell(1, nFile);
for iFile = 1:nFile
  load(fullfile(filePath(iFile).folder, filePath(iFile).name))
  data_cell{iFile} = data;
  params_cell{iFile} = params;
end

data = [data_cell{:}];
params = [params_cell{:}];

figure
subplot(2,1,1)
title('input parameter')
histogram(log10([params(:).a]), 20)
xlabel('Logarithmic infection rate (log(a))')
ylabel('Count')
% set(gca, 'XScale', 'log')
subplot(2,1,2)
title('output response')
histogram([data(:).V_peakTime], 20)
xlabel('Peak time (days)')
ylabel('Count')