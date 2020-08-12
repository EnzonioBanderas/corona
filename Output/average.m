clear all
close all

% Collect all data
folder = '0812_1_1000_sims';
files = dir([folder, '/Gillespie*']);
Nsim = length(files);
alldata = struct();

for vb = 1:Nsim
    fname = [folder,'/', files(vb).name];
    load(fname);
    
    alldata(vb).('time') = titer.data(:,1);
    alldata(vb).('ntot') = titer.data(:,2);
    alldata(vb).('U') = titer.data(:,4);
    alldata(vb).('I') = titer.data(:,5);
    alldata(vb).('V') = titer.data(:,6);
    alldata(vb).('maxD') = maxD;
    alldata(vb).('maxR') = maxR;
    alldata(vb).('statY') = statY;
    alldata(vb).('statR') = statR;
    alldata(vb).('statD') = statD;
end

%% figure 1: time course of infection

orange = [204 102 0];
orange = orange / max(orange);

figure(1)

subplot(2,1,1)
for vb = 1:Nsim
    t = alldata(vb).('time');
    U = alldata(vb).('U');
    V = alldata(vb).('V');
    I = alldata(vb).('I');
    
    semilogy(t, U,'b', 'LineWidth',2)
    hold on
    semilogy(t, V,'r', 'LineWidth',2)
    hold on
    semilogy(t, I,'Color',orange,'LineWidth',2)
    
    hold on
end
semilogy(ti, x(:,1), 'k', 'LineWidth',2)
hold on
semilogy(ti, x(:,2), 'k', 'LineWidth',2)
hold on
semilogy(ti, x(:,3), 'k', 'LineWidth',2)
hold off
legend({'Not infected cells','Free viral particles','Infected cells'},'FontSize',9)
leg.ItemTokenSize = [15, 1];
%title('Viral load vs time')
ylabel('Number of cells or viral propagules')
set(gca,'LineWidth',1)


alltime = vertcat(alldata(:).time)';
allntot = vertcat(alldata(:).ntot)';

[sortedTime, ind] = sort(alltime);
sortedntot = allntot(ind);
Npoints = Nsim;
[timeAxis, ntot_mean, ntot_error] = timeAverage(sortedTime, sortedntot, Npoints);


subplot(2,1,2)
for vb = 1:Nsim
    t = alldata(vb).('time');
    ntot = alldata(vb).('ntot');
    
    plot(t, ntot,'k', 'LineWidth',2)
    hold on
end

plot(timeAxis, ntot_mean, 'r','LineWidth',2)
ylabel('Number of distinct nucleotide sequences')
xlabel('Time (days p.i.)')

set(gca,'LineWidth',1)

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Position',[200 200 900 550],'Units','inches')
saveas(gcf,'Figures/sequences.pdf')

%% Figure 2: mutation analysis]

shadyred = [250,189,189];
shadyred = shadyred / max(shadyred);
red = {shadyred, 'r'};

allmaxR = horzcat(alldata(:).maxR);
sortedMaxR = allmaxR(ind);

Npoints = Nsim;
[timeAxis, maxR_mean, maxR_error] = timeAverage(sortedTime, sortedMaxR, Npoints);

nans = isnan(maxR_mean);
maxR_mean(nans) = [];
maxR_error(nans) = [];
timeAxis(nans) = [];

figure(2)
subplot(2,1,2)

plotShadyError(timeAxis, maxR_mean, maxR_error, red)
ylabel('Max reproduction rate  in population')
xlabel('Time (days p.i.)')
%title('Benificial mutations')
set(gca, 'lineWidth', 1);
set(gca, 'FontSize', 14)

maxY = 0;
for vb = 1:Nsim
    lengthY = length(alldata(vb).statY);
    if maxY < lengthY
        maxY = lengthY;
    end
end
        
statY = zeros(Nsim, maxY);
for vb = 1:Nsim
    if sum(isnan(alldata(vb).statY)) == 0
        statY_vb = alldata(vb).statY;
        statY_vb = [statY_vb, zeros(1, maxY - length(statY_vb))]; % add zeros to make the right length;
        statY(vb, :) = statY_vb;
    else
        disp(vb)
    end
end
% Plot histogram
subplot(2,1,1)
barcolor = [.2 .6 .5];
maxY = find(sum(statY)>0, 1, 'last');
statY = statY(:, 1:maxY);
x = linspace(0, maxY-1, maxY);

meanY = mean(statY);
errorY = std(statY);
bar(x, meanY,'FaceColor', barcolor)
set(gca,'yscale','log')

hold on

er = errorbar(x,meanY, errorY);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

%title('Relative abundance of hamming distances')
xlabel('Hamming distance s')
ylabel('$\bar{N}(s)$','Interpreter','Latex')

set(gca, 'FontSize', 14)

set(gcf,'PaperOrientation','landscape');

set(gcf,'Color','w','Position',[200 200 900 550],'Units','inches')
saveas(gcf,'Figures/mutations.pdf')

%% Figure 3: max values, figure 4: stationary values.
figure()
maxR_all = zeros(1, Nsim); 
maxD_all = zeros(1, Nsim);
statR_all = zeros(1, Nsim);
statD_all = zeros(1, Nsim);

for i = 1:Nsim
    maxR_all(i) = max(alldata(i).maxR);
    maxD_all(i) = max(alldata(i).maxD);
    statR_all(i) = alldata(i).statR;
    statD_all(i) = alldata(i).statD;
end
subplot(2,1,1)
histogram(maxR_all)
xlabel('Max replication rate during the infection')
ylabel('Frequency (out of 1000 simulations)')

subplot(2,1,2)
histogram(maxD_all)
xlabel('Max hamming distance during the infection')
ylabel('Frequency (out of 1000 simulations)')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Position',[200 200 900 550],'Units','inches')
saveas(gcf,'Figures/maxRD_dist.pdf')

figure()
subplot(2,1,1)
histogram(statR_all)
xlabel('Stationary replication rate during the infection')
ylabel('Frequency (out of 1000 simulations)')

subplot(2,1,2)
histogram(statD_all)
xlabel('Stationary hamming distance during the infection')
ylabel('Frequency (out of 1000 simulations)')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Position',[200 200 900 550],'Units','inches')
saveas(gcf,'Figures/statRD_dist.pdf')




function [timeAxis, means, error] = timeAverage(t, y, Npoints)
    time = [t(1:Npoints:end), max(t)];
    timeAxis = (time(2:end) + time(1:end-1))/2;
    
    means = zeros(1,length(timeAxis));
    error = zeros(1,length(timeAxis));

    for i = 1:length(timeAxis)
        inBin = (t > time(i) & t < time(i+1));
        means(i) = nanmean(y(inBin));
        error(i) = nanstd(y(inBin));
    end
end

function plotShadyError(x, means, error, color)

    curve1 = means + error;
    curve2 = max(means - error, 0.00001);
    disp(curve1)
    disp(curve2)
    xaxis = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    
    fill(xaxis, inBetween, color{1});
    hold on;
    plot(x, means, color{2}, 'LineWidth', 2);
    set(gca, 'YScale', 'log' )
end


function [t, x] = solveDE(tRange, initialCond, params)
    
    a = params(1);
    b = params(2);
    c = params(3);
    r0 = params(4);
    
    dfdt = @(t,x) [r0 * x(2) - b * x(1) - a * x(3) * x(1); ... 
                   a * x(3) * x(1) - c * x(2); ... 
                   - a * x(3) * x(1)];
                       
    [t,x] = ode45(dfdt, tRange, initialCond);
end
