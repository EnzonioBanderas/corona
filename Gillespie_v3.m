clear all 
close all

%% Viral evolution with Gillespie algorithm

% Source of model: 
% Woo & Reifman (2013). Quantitative Modeling of Virus Evolutionary
% Dynamics and Adaptation in Serial Passages Using Empirically Inferred
% Fitness Landscapes. Journal of Virology, V88 N2, p1039 - 1050.

% Possible events/reactions:
% Infection: U + Vn -> In     (R1) rate = a , reaction_rate = a*U*Vn 
% Replication: In -> In + Vm  (R2) rate = r(d), reaction_rate = r*I
% Death : In -> 0             (R3) rate = b
% Clearence: Vn -> 0          (R4) rate = b

%% Get protein and genome reference sequences #############################
clear all
[gRefSeq, pRefSeq, pNames, pInfo] = getRefSeq();
[beta0, beta1] = logisticRegressionProteins();

for i = 1:length(pNames)
    protein = pNames{i};
    pInfo.(protein).betas(1) = beta0(i);
    pInfo.(protein).betas(2) = beta1(i);
end
%% Initialize #############################################################
La = length(pRefSeq);                   
L = length(gRefSeq);        % length of genome sequence
a = 2e-5;                   % infection rate per day.
b = 0.4;                    % death/clearence rate of infected cells per day.
V0 = 400;                   % initial number of viruses
U0 = 1e5;                   % initial number of not-infected cells
ksi = 1.0;                 	% fitness decay
sigma = 0.1;                % standard deviation of fitness
T = inf;                    % maximal time (days)

mu_array = 0;            % mutation rate(s). Can be a single float (one mutation rate) or an array.
mu = 1e-3;

%WTseq = 'catacacaagga';        % initial genome sequence
%refAseq = nt2aa(WTseq);        % reference aa sequence (with highest fitness). 
d0 = zeros(length(pNames),1);  	% distance of WTseq to fittest strain
r0 = 2;                         % Fitness of reference sequence 

infection = boolean(zeros(1,length(mu_array)));     % store whether a full infection occurred
Uarray = zeros(1,length(mu_array));                 % store the number of not-infected cells at the end of the infection
titer = struct();
Nrep = 1;

% Start for-loop over mutation rates ######################################
for k = 1:Nrep
    % initialize
    disp(k)
    U = U0;                                                                 % number of infected cells.
    seq = {}; aseq = {};
    seq{1} = gRefSeq;                                                     	% cell array to keep track of nt sequences
    aseq{1} = pRefSeq;                                                      % cell array to keep track of aa sequences   
    % counters
    ntot = 1;                                                               % initial number of distinct viral genotypes
    nAA = 1;                                                                % initial number of distinct viral phenotypes
    s = 1;                                                                  % counter for while loop
    m = 1;                                                                  % counter to collect stats
    t = 0.0;                                                                % initial time (days)
    % arrays to collect stats
    V = [V0];                                                               % array for number of viruses of all strains
    I = [0];                                                                % array for number of infected cells of all strains
    Y = [V(1), zeros(1, L)];                                                % matrix for number of viruses with distance d from WT seq.
    r = [r0];                                                               % array for fitnesses 
    d = [d0];                                                               % array for distances
    data = [t, ntot, nAA, U, sum(I), sum(V)];                             	% matrix to collect time, # distinct genotypes, # distinct phenotypes, # not-infected cells, # infected cells, # viruses.
    meanDistance = [0];                                                    	% mean distances
    meanFitness = [r0];                                                   	% mean fitness
        
    % Start while loop ####################################################
    while t < T
        if sum(I+V) == 0                                                  	% Stop if there are no more viral particles left
             break
        end            
         s = s + 1;                                                         % counter
         % calculate total reaction rate
         Rtot = (a*U + b)*sum(V) + b*sum(I) + sum(r.*I);                                          
         % take time step
         dt = (1/Rtot) * log(1/rand);                                       % time step is exponentially distributed and depends on total reaction rate
         t = t + dt;
         % choose a random number between 0 and Rtot
         rnd = Rtot * rand;
         c = 0;
         % loop over viral strains
         for i = 1:ntot    
             % choose which of the 4 reactions will occur 
             % and for which strain
             c = c + a*U*V(i);                                              % R1: infection by strain i
             if c > rnd                                                         
                 I(i) = I(i) + 1;                                               
                 V(i) = V(i) - 1;
                 U = U - 1;                                                    
                 break
             end
             c = c + r(i)*I(i);                                             % R2: replication of strain i
             if c > rnd
                 i0 = i;
                 [seq, aseq, ntot, V, r, I, d] = ...
                     replicate(i, seq, aseq, mu, ntot, V, r, I, d, ... 
                     sigma,r0, pRefSeq, pInfo);
                 break
             end
             c = c + b*I(i);                                                % R3: clearence of a cell infected by strain i
             if c > rnd
                 I(i) = I(i) - 1;
                 break
             end
             c = c + b*V(i);                                                % R4: clearence of a free viral particle of strain i
             if c > rnd
                 V(i) = V(i) - 1;
                 break
             end
         end % for-loop
         % remove strains that have gone extinct
         [seq, aseq, V, I, r, d, ntot] = ...
             remove(seq, aseq, V, I, r, d, ntot);
         % collect statistics every 100 iterations and at the end
         if mod(s,100) == 0
             nAA = length(unique(r));
             m = m + 1;             
             for j = 1:size(Y,2)
                dist = j - 1;
                d_logical = (sum(d,1) == dist);
                sumV = sum(V(d_logical));
                sumI = sum(I(d_logical));
                Y(m,j) = sumV + sumI;
             end
             meanFitness(m) = sum(r.*(V+I)) / sum(V+I);
             meanDistance(m) = sum(sum(d,1).*(V+I)) / sum(V+I);
             data(m,:) = [t, ntot, nAA, U, sum(I), sum(V)];
         end
         % display progression every 10000 iterations
         if mod(s,100) == 0
             disp(['t=',num2str(t),', ntot=',num2str(ntot)])
         end
    end % while-loop
    titer(k).time = data(:,1);
    titer(k).load = data(:,5) + data(:,6);
    % Has infection occured?
    if U == 0
        infection(k) = true;
    end
    % Update error threshold detectors:
    time = data(:,1)';
    delta_t = time(2:end) - time(1:end-1);
    sumY = 1 ./ sum(Y,2);
    relativeY = Y .* repmat(sumY, 1, size(Y,2));
    statY(k,:) = (delta_t * relativeY(2:end, :)) / time(end) ;
    statD(k) = sum(meanDistance(2:end) .* delta_t) / time(end);
    statR(k) = sum(meanFitness(2:end) .* delta_t) / time(end); 
    Uarray(k) = U;
end

save('Gillespie.mat','-v7.3');

%% plot the results #######################################################
% if length(mu_array) == 1    % if only one mutation rate:
%     
%     if Nrep == 1
%     
%         figure()
%         subplot(2,1,1)
%         semilogy(dataP(:,1), dataP(:,4), dataP(:,1), dataP(:,5), dataP(:,1), dataP(:,6),'LineWidth',2)
%         xlabel('time (days)')
%         ylabel('Cell count')
%         legend('Number of uninfected cells', 'Number of infected cells','Number of free viral particles')
%         title('Viral load vs time')
%         set(gca, 'Fontsize', 24, 'LineWidth',2)
% 
%         subplot(2,1,2)
%         plot(dataP(:,1), dataP(:,2), dataP(:,1), dataP(:,3),'LineWidth',2)
%         xlabel('Time (days)')
%         ylabel('Number of viral strains')
%         title('Number of viral strains vs time')
%         legend('Number of distinct nucleotide sequences', 'Number of distinct AA sequences')
%         set(gca, 'Fontsize', 24, 'LineWidth',2)
% 
%         figure()
%         subplot(2,1,1)
%         plot(dataP(:,1), Y)
%         xlabel('time (days)')
%         ylabel('Viral load')
%         presentDistances = find(sum(Y,1)~=0);
%         for i = 1:length(presentDistances)
%             leg{i} = ['d=',num2str(presentDistances(i)-1)];
%         end
%         legend(leg)
%         title('Abundance of strains with distance d from the WT sequence')
%         set(gca, 'Fontsize', 24)
% 
%         subplot(2,1,2)
%         plot(dataP(:,1), Y(:,1), dataP(:,1), sum(Y(:,2:end),2))
%         legend('d=0','error tail')
%         xlabel('time (days)')
%         ylabel('Viral load')
%         title('Original sequences versus error tail')
%         set(gca, 'Fontsize', 24)
%         
%     else
%         figure()
%         for i=1:length(titer)
%             plot(titer(i).time, titer(i).load)
%             hold on
%         end
%         hold off
%     end
%     
% else                    % if mutiple mutation rates:
%     figure()
%     dist0 = statY(:,1);
%     dist1 = statY(:,2);
%     dist2 = statY(:,3);
%     dist3 = statY(:,4);
%     error_tail = sum(statY(:,5:end),2);
%     subplot(2,2,1:2)
%     plot(mu_array, [dist0, dist1, dist2, dist3, error_tail], mu_array,(U0-Uarray), '--')
%     leg = {'d=0','d=1','d=2','d=3','d>3','Number of infected cells'};
%     legend(leg,'Location','East','Fontsize',16)
%     title('Stationary value of particle number','Fontsize',16)
%     xlabel('Mutation rate \mu','Fontsize',16,'Fontsize',16)
%     ylabel('$\bar{N_d}(\mu)$','Interpreter','latex','Fontsize',16)
%     text(0.1,0.8,{'\leftarrow infection occurs (=1)','or not (=0)'})
%     text(0.2,0.4,'$\bar{N_d}(\mu) = \frac{1}{T}\int_{t=0}^TN_d(t,\mu)dt$',...
%         'Interpreter','latex','Fontsize',16)
% 
%     subplot(2,2,3)
%     plot(mu_array, statD,'o')
%     title('Stationary value of Hamming distance','Fontsize',16)
%     xlabel('Mutation rate \mu','Fontsize',16)
%     ylabel('$\bar{d}(\mu)$','Interpreter','latex','Fontsize',16)
% 
%     subplot(2,2,4)
%     plot(mu_array, statR,'o')
%     title('Stationary value of mean fitness','Fontsize',16)
%     xlabel('Mutation rate \mu','Fontsize',16)
%     ylabel('$\bar{R}(\mu)$','Interpreter','latex','Fontsize',16)
% 
%     suptitle('Error catastrophe with Gillespie algorithm')
% end
