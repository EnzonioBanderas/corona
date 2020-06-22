clear all
% close all
% clc

load('Gillespie.mat');

% data [t, ntot, nAA, U, sum(I), sum(V)]
data = vertcat(titer(:).data); % before doing this individually average values within time bins! In that way SEM is determined by iteration number, not by sampling resolution

% group time intervals
timeStep = 0.25;
timeInt = round(data(:,1)/timeStep);
time = unique(timeInt)*dt;
timeGroups = findgroups(timeInt);
% [timeGroups, sortVec] = sort(timeGroups);
% data = data(sortVec, :);

% number of genotypic variants
ntot_mean = splitapply(@mean, data(:,2), timeGroups);
ntot_std = splitapply(@std, data(:,2), timeGroups);
ntot_length = splitapply(@length, data(:,2), timeGroups);
ntot_sem = ntot_std ./ sqrt(ntot_length);

% number of uninfected cells
U_mean = splitapply(@mean, data(:,4), timeGroups);
U_std = splitapply(@std, data(:,4), timeGroups);
U_length = splitapply(@length, data(:,4), timeGroups);
U_sem = U_std ./ sqrt(U_length);

% sum of infected cells
I_sum_mean = splitapply(@mean, data(:,5), timeGroups);
I_sum_std = splitapply(@std, data(:,5), timeGroups);
I_sum_length = splitapply(@length, data(:,5), timeGroups);
I_sum_sem = I_sum_std ./ sqrt(I_sum_length);

% sum of viral particles
V_sum_mean = splitapply(@mean, data(:,6), timeGroups);
V_sum_std = splitapply(@std, data(:,6), timeGroups);
V_sum_length = splitapply(@length, data(:,6), timeGroups);
V_sum_sem = V_sum_std ./ sqrt(V_sum_length);



figure()
subplot(2,1,1)
c = colormap('lines');
hold on
plot(time, U_mean, 'Color', c(1,:), 'linewidth', 2)
plot(time, I_sum_mean, 'Color', c(2,:), 'linewidth', 2)
plot(time, V_sum_mean, 'Color', c(3,:), 'linewidth', 2)
% figure
% V_sum_std_logical = V_sum_std~=0;
% V_sum_std_plot = V_sum_std(V_sum_std_logical);
% time_plot = time(V_sum_std_logical);
legend('Number of uninfected cells', 'Number of infected cells','Number of free viral particles','AutoUpdate','off')
patch_y = [U_mean + U_sem; ...
           U_mean(end) - U_sem(end); ...
           U_mean(end:-1:1) - U_sem(end:-1:1); ...
           U_mean(1) + U_sem(1)];
patch_y(patch_y == 0) = 0.01;
patch([time; time(end); time(end:-1:1); time(1)], ...
     patch_y, ...
     c(1,:), 'FaceAlpha', 0.3)
patch_y = [I_sum_mean + I_sum_sem; ...
           I_sum_mean(end) - I_sum_sem(end); ...
           I_sum_mean(end:-1:1) - I_sum_sem(end:-1:1); ...
           I_sum_mean(1) + I_sum_sem(1)];
patch_y(patch_y == 0) = 0.01;
patch([time; time(end); time(end:-1:1); time(1)], ...
     patch_y, ...
     c(2,:), 'FaceAlpha', 0.3)
 patch_y = [V_sum_mean + V_sum_sem; ...
           V_sum_mean(end) - V_sum_sem(end); ...
           V_sum_mean(end:-1:1) - V_sum_sem(end:-1:1); ...
           V_sum_mean(1) + V_sum_sem(1)];
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
plot(time, ntot_mean, 'Color', c(1,:), 'linewidth', 2)
 patch_y = [ntot_mean + ntot_sem; ...
           ntot_mean(end) - ntot_sem(end); ...
           ntot_mean(end:-1:1) - ntot_sem(end:-1:1); ...
           ntot_mean(1) + ntot_sem(1)];
patch_y(patch_y == 0) = 0.01;
patch([time; time(end); time(end:-1:1); time(1)], ...
     patch_y, ...
     c(1,:), 'FaceAlpha', 0.3)
xlabel('Time (days)')
ylabel('Number of viral strains')
title('Number of viral strains vs time')
legend('Number of distinct nucleotide sequences')
% legend('Number of distinct nucleotide sequences', 'Number of distinct AA sequences')
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