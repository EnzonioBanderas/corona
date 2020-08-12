clear all
close all
clc

load('Gillespie.mat');
t = 'R = prod(ri)';
mu_array = 1e-6;
data = titer(1).data;
%%
if length(mu_array) == 1    % if only one mutation rate:
    
    figure()
    subplot(2,1,1)
    semilogy(data(:,1), data(:,4), data(:,1), data(:,5), data(:,1), data(:,6), 'LineWidth',2)
    xlabel('time (days)')
    ylabel('Cell count')
    L = {'Number of uninfected cells', 'Number of infected cells','Number of free viral particles'};
    leg = legend(L, 'Fontsize',9);
    leg.ItemTokenSize = [15, 1];
    title('Viral load vs time')
    set(gca,'LineWidth',1)
    
   	subplot(2,1,2)
    plot(data(:,1), data(:,2), data(:,1), data(:,3),'LineWidth',2)
    xlabel('Time (days)')
    ylabel('Number of viral strains')
    title('Number of viral strains vs time')
    L = {'Number of distinct nucleotide sequences', 'Number of distinct AA sequences'};
    leg = legend(L, 'Fontsize',9);
    leg.ItemTokenSize = [15, 1];
    set(gca,'LineWidth',1)
    
    suptitle(t)
    set(gcf,'PaperOrientation','landscape');
    saveas(gcf,'Figures/load.pdf')

    figure()
    subplot(2,1,1)
    plot(data(:,1), Y(:,1:20), 'LineWidth',2)
    xlabel('time (days)')
    ylabel('Viral load')
    presentDistances = find(sum(Y,1)~=0);
%    L = {};
%     for i = 1:length(presentDistances)
%         L{i} = ['d=',num2str(presentDistances(i)-1)];
%     end
    %legend(L)
    title('Abundance of strains with distance d from the WT sequence')
    set(gca,'LineWidth',1)

    subplot(2,1,2)
    plot(data(:,1), Y(:,1), data(:,1), sum(Y(:,2:end),2), 'LineWidth',2)
    L = {'d=0','error tail'};
    leg = legend(L, 'Fontsize',9);
    leg.ItemTokenSize = [15, 1];
    xlabel('time (days)')
    ylabel('Viral load')
    title('Original sequences versus error tail')
    set(gca,'LineWidth',1)
    
    suptitle(t)
    set(gcf,'PaperOrientation','landscape');
    saveas(gcf,'Figures/strains.pdf')
    
    figure()
    
    subplot(2,1,1)
    plot(data(:,1), meanDistance,'LineWidth',2)
    xlabel('Time (days)')
    ylabel('Mean distance')
    title('Weighted mean distance of the population')
    set(gca, 'Linewidth',1)
    set(gcf,'PaperOrientation','landscape');
    
    subplot(2,1,2)
    plot(data(:,1), meanFitness,'LineWidth',2)
    xlabel('Time (days)')
    ylabel('Mean fitness')
    title('Weighted mean reproduction rate of the population')
    set(gca, 'Linewidth',1)
    set(gcf,'PaperOrientation','landscape');
    
    suptitle(t)
    set(gcf,'Color','w','Position',[250 250 800 400],'Units','inches')
    saveas(gcf,'Figures/meanReprRate_Distance.pdf')
    
    figure()
    subplot(2,1,1)
    plot(data(:,1), maxD,'LineWidth',2)
    xlabel('Time (days)')
    ylabel('Max distance')
    title('Max distance from wildtype sequence')
    set(gca, 'Linewidth',1)
    
    subplot(2,1,2)
    plot(data(:,1), maxR,'LineWidth',2)
    xlabel('Time (days)')
    ylabel('Max reproduction rate')
    title('Max reproduction rate in population')
    set(gca, 'Linewidth',1)
   
   	suptitle(t)
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Color','w','Position',[250 250 800 400],'Units','inches')
    saveas(gcf,'Figures/maxReprRate_Distance.pdf')

    figure()
    maxM = size(Y,1);
    l_array = round(linspace(maxM/6, maxM, 6));
    for n = 1:length(l_array)
        l = l_array(n);
        histDist = [];
        for i=1:find(Y(l,:),1,'last')
            histDist = [histDist, repmat(i-1, [1, Y(l,i)])];
        end
        subplot(2,3,n)
        histogram(histDist,'Normalization','pdf');
        ylim([0,1])
        xlim([0,max(maxD)])
        xlabel('Number of mutations')
        ylabel('Relative abundance')
        title(['t = ',num2str(data(l,1)),' days'])
    end
    suptitle([t,' histograms of mutation number'])
    set(gcf,'PaperOrientation','landscape');
	set(gcf,'Color','w','Position',[250 250 800 400],'Units','inches')
    saveas(gcf,'Figures/mutationHist.pdf')

    
    
else                    % if mutiple mutation rates:
    figure()
    dist0 = statY(:,1);
    dist1 = statY(:,2);
    error_tail = sum(statY(:,3:end),2);
    subplot(2,2,1:2)
    plot(mu_array, [dist0, dist1, error_tail], mu_array, (U0-Uarray)/U0, '--', 'LineWidth',2)
    xlabel('Mutation rate \mu','Fontsize',12)
    ylabel('Rel. particle number','Fontsize',12)
    L = {'d=0','d=1','d>1','Infected cells'};
    
    leg = legend(L,'Location','East','Fontsize',9);
    leg.ItemTokenSize = [15, 1];
    
    set(gca, 'Linewidth',1)

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

    suptitle(t)
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Color','w','Position',[250 250 800 500],'Units','inches')
    saveas(gcf,'Figures/mutationRates.pdf')
end

%%


% hold on
% pd = fitdist(histDist'+1, 'GeneralizedPareto');
% x = linspace(0,6,100);
% pdf_ = pdf(pd, x);
% plot(x, pdf_)
