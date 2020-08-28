clear all



% N = 100;
% for x = 1:N
%     x_str = num2str(x);
%     Gillespie_func(x_str);
% end



filePath = dir('a_linear_uncertainty/a_linear_uncertainty*');
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

[input, output_sort] = sort([params(:).a]);

figure
subplot(2,1,1)
title('input parameter')
histogram(log10([params(:).a]), 20)
xlabel('Logarithmic infection rate (log(a))')
ylabel('Count')
% set(gca, 'XScale', 'log')
subplot(2,1,2)
title('output response')
histogram([data(output_sort).V_peakTime], 20)
xlabel('Peak time (days)')
ylabel('Count')

figure

output_cell = {[data(output_sort).V_peakTime], ...
               [data(output_sort).V_peak], ...
               [data(output_sort).t_end]};
output_str = {['V peak time (', num2str(max(output_cell{1})), ')'], ...
              ['V peak (', num2str(max(output_cell{2})), ')']', ...
              ['t end (', num2str(max(output_cell{3})), ')']'};

for i=1:length(output_cell)
    plot([params(output_sort).a], output_cell{i}/max(output_cell{i})), hold on
end

set(gca, 'XScale', 'log')
legend(output_str)
xlabel('Infection rate')
ylabel('Output max ratio')






figure

output_cell = {[data(output_sort).V_peakTime], ...
               [data(output_sort).t_end]};
output_str = {'V peak time', ...
              't end'};

for i=1:length(output_cell)
    plot([params(output_sort).a], output_cell{i}), hold on
end

set(gca, 'XScale', 'log')
legend(output_str)
xlabel('Infection rate')
ylabel('Output max ratio')