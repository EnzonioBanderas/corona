clear all
% close all
% clc

load('Gillespie.mat');

% data [t, ntot, nAA, U, sum(I), sum(V)]
timeStep = 0.25;
nIter = length(titer);
nMeasure = size(titer(1).data, 2);
for k = 1:nIter
    data = titer(k).data;
    
    timeInt = round(data(:,1)/timeStep);
    time = unique(timeInt)*timeStep;
    timeGroups = findgroups(timeInt);
    
    % Fill in data_mean values
    data_mean = zeros(length(time), nMeasure);
    data_mean(:, 1) = time;
    for iMeasure = 2:nMeasure
        data_mean(:, iMeasure) = splitapply(@mean, data(:,iMeasure), timeGroups);
%         ntot_std = splitapply(@std, data(:,iMeasure), timeGroups);
%         ntot_length = splitapply(@length, data(:,iMeasure), timeGroups);
%         data_sem(:, 2) = ntot_std ./ sqrt(ntot_length);
    end
    
    % Assign data_mean to iteration structure
    titer(k).data_mean = data_mean;
end
data = vertcat(titer(:).data_mean); % concatenate meaned data

% Repeat process but also calculate SEM
timeInt = round(data(:,1)/timeStep);
time = unique(timeInt)*timeStep;
timeGroups = findgroups(timeInt);

% Fill in data_mean, data_std, data_length and data_sem values
data_mean = zeros(length(time), nMeasure); data_mean(:, 1) = time;
data_std = zeros(length(time), nMeasure); data_std(:, 1) = time;
data_length = zeros(length(time), nMeasure); data_length(:, 1) = time;
data_sem = zeros(length(time), nMeasure); data_sem(:, 1) = time;
for iMeasure = 2:nMeasure
    data_mean(:, iMeasure) = splitapply(@mean, data(:,iMeasure), timeGroups);
    data_std(:, iMeasure) = splitapply(@std, data(:,iMeasure), timeGroups);
    data_length(:, iMeasure) = splitapply(@length, data(:,iMeasure), timeGroups);
    data_sem(:, iMeasure) = data_std(:, iMeasure) ./ sqrt(data_length(:, iMeasure));
end



figure()
subplot(2,1,1)
c = colormap('lines');
hold on
plot(time, data_mean(:, 4), 'Color', c(1,:), 'linewidth', 2)
plot(time, data_mean(:, 5), 'Color', c(2,:), 'linewidth', 2)
plot(time, data_mean(:, 6), 'Color', c(3,:), 'linewidth', 2)
% figure
% V_sum_std_logical = V_sum_std~=0;
% V_sum_std_plot = V_sum_std(V_sum_std_logical);
% time_plot = time(V_sum_std_logical);
legend('Number of uninfected cells', 'Number of infected cells','Number of free viral particles','AutoUpdate','off')
patch_y = [data_mean(:, 4) + data_sem(:, 4); ...
           data_mean(end, 4) - data_sem(end, 4); ...
           data_mean(end:-1:1, 4) - data_sem(end:-1:1, 4); ...
           data_mean(1, 4) + data_sem(1, 4)];
patch_y(patch_y == 0) = 0.01;
patch([time; time(end); time(end:-1:1); time(1)], ...
     patch_y, ...
     c(1,:), 'FaceAlpha', 0.3)
patch_y = [data_mean(:, 5) + data_sem(:, 5); ...
           data_mean(end, 5) - data_sem(end, 5); ...
           data_mean(end:-1:1, 5) - data_sem(end:-1:1, 5); ...
           data_mean(1, 5) + data_sem(1, 5)];
patch_y(patch_y == 0) = 0.01;
patch([time; time(end); time(end:-1:1); time(1)], ...
     patch_y, ...
     c(2,:), 'FaceAlpha', 0.3)
 patch_y = [data_mean(:, 6) + data_sem(:, 6); ...
           data_mean(end, 6) - data_sem(end, 6); ...
           data_mean(end:-1:1, 6) - data_sem(end:-1:1, 6); ...
           data_mean(1, 6) + data_sem(1, 6)];
patch_y(patch_y == 0) = 0.01;
patch([time; time(end); time(end:-1:1); time(1)], ...
     patch_y, ...
     c(3,:), 'FaceAlpha', 0.3)
xlabel('time (days)')
ylabel('Cell count')
title('Viral load vs time')
set(gca, 'YScale', 'log')
set(gca, 'Fontsize', 24)

subplot(2,1,2)
hold on
plot(time, data_mean(:, 2), 'Color', c(1,:), 'linewidth', 2)
plot(time, data_mean(:, 3), 'Color', c(2,:), 'linewidth', 2)
legend('Number of distinct nucleotide sequences', 'Number of distinct AA sequences', 'AutoUpdate','off')
patch_y = [data_mean(:, 2) + data_sem(:, 2); ...
           data_mean(end, 2) - data_sem(end, 2); ...
           data_mean(end:-1:1, 2) - data_sem(end:-1:1, 2); ...
           data_mean(1, 2) + data_sem(1, 2)];
% patch_y(patch_y == 0) = 0.01;
patch([time; time(end); time(end:-1:1); time(1)], ...
     patch_y, ...
     c(1,:), 'FaceAlpha', 0.3)
patch_y = [data_mean(:, 3) + data_sem(:, 3); ...
           data_mean(end, 3) - data_sem(end, 3); ...
           data_mean(end:-1:1, 3) - data_sem(end:-1:1, 3); ...
           data_mean(1, 3) + data_sem(1, 3)];
% patch_y(patch_y == 0) = 0.01;
patch([time; time(end); time(end:-1:1); time(1)], ...
     patch_y, ...
     c(2,:), 'FaceAlpha', 0.3)
xlabel('Time (days)')
ylabel('Number of viral strains')
title('Number of viral strains vs time')
set(gca, 'Fontsize', 24)

% figure()
% subplot(2,1,1)
% plot(time, Y)
% xlabel('time (days)')
% ylabel('Viral load')
% presentDistances = find(sum(Y,1)~=0);
% for i = 1:length(presentDistances)
%     leg{i} = ['d=',num2str(presentDistances(i)-1)];
% end
% %legend(leg)
% title('Abundance of strains with distance d from the WT sequence')
% set(gca, 'Fontsize', 24)
% 
% subplot(2,1,2)
% plot(data(:,1), Y(:,1), data(:,1), sum(Y(:,2:end),2))
% legend('d=0','error tail')
% xlabel('time (days)')
% ylabel('Viral load')
% title('Original sequences versus error tail')
% set(gca, 'Fontsize', 24)