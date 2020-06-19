clear all
% close all
% clc

load('Gillespie.mat');

if length(mu_array) == 1    % if only one mutation rate:
    
    figure()
    subplot(2,1,1)
    semilogy(data(:,1), data(:,4), data(:,1), data(:,5), data(:,1), data(:,6))
    xlabel('time (days)')
    ylabel('Cell count')
    legend('Number of uninfected cells', 'Number of infected cells','Number of free viral particles')
    title('Viral load vs time')
    set(gca, 'Fontsize', 24)
    
   	subplot(2,1,2)
    plot(data(:,1), data(:,2), data(:,1), data(:,3))
    xlabel('Time (days)')
    ylabel('Number of viral strains')
    title('Number of viral strains vs time')
    legend('Number of distinct nucleotide sequences', 'Number of distinct AA sequences')
    set(gca, 'Fontsize', 24)

    figure()
    subplot(2,1,1)
    plot(data(:,1), Y)
    xlabel('time (days)')
    ylabel('Viral load')
    presentDistances = find(sum(Y,1)~=0);
    for i = 1:length(presentDistances)
        leg{i} = ['d=',num2str(presentDistances(i)-1)];
    end
    %legend(leg)
    title('Abundance of strains with distance d from the WT sequence')
    set(gca, 'Fontsize', 24)

    subplot(2,1,2)
    plot(data(:,1), Y(:,1), data(:,1), sum(Y(:,2:end),2))
    legend('d=0','error tail')
    xlabel('time (days)')
    ylabel('Viral load')
    title('Original sequences versus error tail')
    set(gca, 'Fontsize', 24)
    
else                    % if mutiple mutation rates:
    figure()
    dist0 = statY(:,1);
    dist1 = statY(:,2);
    error_tail = sum(statY(:,3:end),2);
    subplot(2,2,1:2)
    plot(mu_array, [dist0, dist1, error_tail], mu_array, (U0-Uarray)/U0, '--', 'LineWidth',2)
    leg = {'d=0','d=1','d>1','Number of infected cells'};
    legend(leg,'Location','East','Fontsize',20)
    title('Average number of viruses','Fontsize',24)
    xlabel('Mutation rate \mu','Fontsize',24)
    ylabel('Avg particle number','Fontsize',24)
    %text(0.2,0.8,'$\bar{N_d}(\mu) = \frac{1}{T}\int_{t=0}^TN_d(t,\mu)dt$',...
    %    'Interpreter','latex','Fontsize',16)
    set(gca, 'Linewidth',2)

    subplot(2,2,3)
    plot(mu_array, statD,'o')
    title('Stationary value of Hamming distance','Fontsize',16)
    xlabel('Mutation rate \mu','Fontsize',16)
    ylabel('$\bar{d}(\mu)$','Interpreter','latex','Fontsize',16)

    subplot(2,2,4)
    plot(mu_array, statR,'o')
    title('Stationary value of mean fitness','Fontsize',16)
    xlabel('Mutation rate \mu','Fontsize',16)
    ylabel('$\bar{R}(\mu)$','Interpreter','latex','Fontsize',16)

    suptitle('Error catastrophe with Gillespie algorithm')
end