clear all 
close all

% Implementation of SIR model and viral evolution model with Gillespie
% algorithm.
% TOTAL RUNTIME ~ 10 mins.

%% SIR model with Gillespie algortithm
% Source: Wikipedia https://en.wikipedia.org/wiki/Gillespie_algorithm

% Initialize
N = 1000;                                                                   % population size
T = 50;                                                                     % Maximum elapsed time
t = 0;                                                                      % start time
V = 100;                                                                    % Spatial parameter 
alpha = 0.5;                                                                % rate of infection after contact
beta = 0.5;                                                                 % rate of cure
n_I = 1;                                                                    % initial infected population

% Compute susceptible popultion, set recovered to zero
n_S = N - n_I;
n_R = 0;

% Start time loop
s = 0;
while t < T
    s = s + 1;
    if n_I == 0                                                             % Stop loop if there are no more individuals infected
        break
    end
    % compute reaction rates
    w1 = alpha * n_S * n_I / V;                                             % first reaction: rate with which a susceptible person gets infected
    w2 = beta * n_I;                                                        % second reaction: rate with which an infected person is cured
    W = w1 + w2;                                                            % total rate
    % compute time step
    dt = -log(rand) / W;                                                    % time step is exponentially distributed and depends on total reaction rate
    t = t + dt;
    % choose which reaction occurs in this time step
    if rand < w1 / W                                                        % first reaction: S -> I                                                  
        n_S = n_S - 1;
        n_I = n_I + 1;
    else                                                                    % second reaction: I -> R
        n_I = n_I - 1;
        n_R = n_R + 1;
    end
    % collect statistics
    SIR_data(s,:) = [t, n_S, n_I, n_R];                                    
end

% Plot 
t = SIR_data(:,1);
n_S = SIR_data(:,2);
n_I = SIR_data(:,3);
n_R = SIR_data(:,4);

figure()
plot(t,n_S,'r', t,n_I,'b', t,n_R,'g')
title('Stochastic SIR model with Gillespie algorithm')
xlabel('Time (days)')
ylabel('Number of individuals')
legend('Susceptible', 'Infected', 'Removed')

%% Viral evolution with Gillespie algorithm

% Source of model: 
% Woo & Reifman (2013). Quantitative Modeling of Virus Evolutionary
% Dynamics and Adaptation in Serial Passages Using Empirically Inferred
% Fitness Landscapes. Journal of Virology, V88 N2, p1039 - 1050.

% Possible events/reactions:
% Infection: U + Vn -> In     (R1)
% Replication: In -> In + Vm  (R2)
% Death : In -> 0             (R3)
% Clearence: Vn -> 0          (R4)

% Initialize ##############################################################
clear all

La = 3.0;                   % length of AA sequence
a = 1.0e-3;                 % infection rate per day.
b = 1.5;                    % death/clearence rate of infected cells per day.
V0 = 400;                   % initial number of viruses
U0 = 1e6;                   % initial number of not-infected cells
ksi = 1.0;                 	% fitness decay
sigma = 0.1;                % standard deviation of fitness
T = 15;                     % maximal time (days)

mu_array = 1e-5;          	% mutation rate(s). Can be a single float (one mutation rate) or an array.

WTseq = 'catacacaagga';     % initial genome sequence
refAseq = 'HTQG';           % reference aa sequence (with highest fitness). 
d0 = 0;                     % distance of WTseq to fittest strain
r0 = 6.0;                   % Fitness of reference sequence 

infection = boolean(zeros(1,length(mu_array)));     % store whether a full infection occurred

% Start for-loop over error rates #########################################
for k = 1:length(mu_array)
    % initialize
    mu = mu_array(k);  
    disp(['mu = ',num2str(mu)]);
    U = U0;                                                                 % number of infected cells.
    seq = {}; aseq = {};
    seq{1} = WTseq;                                                         % cell array to keep track of nt sequences
    aseq{1} = nt2aa(seq{1});                                                % cell array to keep track of aa sequences   
    % counters
    ntot = 1;                                                               % initial number of distinct viral genotypes
    nAA = 1;                                                                % initial number of distinct viral phenotypes
    s = 1;                                                                  % counter for while loop
    m = 1;                                                                  % counter to collect stats
    t = 0.0;                                                                % initial time (days)
    % arrays to collect stats
    V = [V0];                                                               % array for number of viruses of all strains
    I = [0];                                                                % array for number of infected cells of all strains
    Y = [V(1), zeros(1, length(WTseq))];                                    % matrix for number of viruses with distance d from WT seq.
    r = [r0];                                                               % array for fitnesses 
    d = [d0];                                                               % array for distances
    data = [t, ntot, nAA, U, sum(I), sum(V)];                               % matrix to collect time, # distinct genotypes, # distinct phenotypes, # not-infected cells, # infected cells, # viruses.
    meanDistance = [0];                                                    	% mean distances
    meanFitness = [r0];                                                   	% mean fitness
    
    % Start loop ##########################################################
    while t < T
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
                     replicate(i, seq, aseq, mu, ntot, V, r, I, d, ksi, sigma,r0, refAseq);
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
         if mod(s,100) == 0 || sum(I+V) == 0
             nAA = length(unique(aseq));
             m = m + 1;             
             for j = 1:size(Y,2)
                dist = j - 1;
                sumV = sum(V(d == dist));
                sumI = sum(I(d == dist));
                Y(m,j) = sumV + sumI;
             end
             meanFitness(m) = sum(r.*(V+I)) / sum(V+I);
             meanDistance(m) = sum(d.*(V+I)) / sum(V+I);
             data(m,:) = [t, ntot, nAA, U, sum(I), sum(V)];
         end
         % display progression every 10000 iterations
         if mod(s,1000) == 0
             disp(['t=',num2str(t),', ntot=',num2str(ntot)])
         end
         if sum(I+V) == 0                                                  	% Stop if there are no more viral particles left
             break
         end
    end % while-loop
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
end

% plot the results #######################################################
if length(mu_array) == 1    % if only one mutation rate:
    
    figure()
    subplot(2,1,1)
    semilogy(data(:,1), data(:,4), data(:,1), data(:,5), data(:,1), data(:,6))
    xlabel('time (days)')
    ylabel('Cell count')
    legend('Number of uninfected cells', 'Number of infected cells','Number of free viral particles')
    title('Viral load vs time')
    
   	subplot(2,1,2)
    semilogy(data(:,1), data(:,2), data(:,1), data(:,3))
    xlabel('Time (days)')
    ylabel('Number of viral strains')
    title('Number of viral strains vs time')
    legend('Number of distinct nucleotide sequences', 'Number of distinct AA sequences')

    figure()
    subplot(2,1,1)
    plot(data(:,1), Y)
    xlabel('time (days)')
    ylabel('Viral load')
    for i = find(sum(Y,1)~=0)
        leg{i} = ['d=',num2str(i-1)];
    end
    legend(leg)
    title('Abundance of strains with distance d from the WT sequence')

    subplot(2,1,2)
    plot(data(:,1), Y(:,1), data(:,1), sum(Y(:,2:end),2))
    legend('d=0','error tail')
    xlabel('time (days)')
    ylabel('Viral load')
    title('Original sequences versus error tail')
    
else                    % if mutiple mutation rates:
    figure()
    dist0 = statY(:,1);
    dist1 = statY(:,2);
    dist2 = statY(:,3);
    dist3 = statY(:,4);
    error_tail = sum(statY(:,5:end),2);
    subplot(2,2,1:2)
    plot(mu_array, [dist0, dist1, dist2, dist3, error_tail], mu_array, infection, '--')
    leg = {'d=0','d=1','d=2','d=3','d>3'};
    legend(leg,'Location','East','Fontsize',16)
    title('Stationary value of particle number','Fontsize',16)
    xlabel('Mutation rate \mu','Fontsize',16,'Fontsize',16)
    ylabel('$\bar{N_d}(\mu)$','Interpreter','latex','Fontsize',16)
    text(0.1,0.8,{'\leftarrow infection occurs (=1)','or not (=0)'})
    text(0.2,0.4,'$\bar{N_d}(\mu) = \frac{1}{T}\int_{t=0}^TN_d(t,\mu)dt$',...
        'Interpreter','latex','Fontsize',16)

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

%% Functions 
function r = fitness(d,r0, ksi, sigma)
% This function determines the fitness of a strain with a distance d from
% the reference sequence, which has fitness r0.
% The fitness is normal distributed with mean mu = r0*exp(-d/ksi) and
% standard deviation sigma.

    %Normal distributed fitness:
    mu = r0*exp(-d/ksi);        
    r = sigma*randn + mu;

    % Single-peak fitness:
%     if d==0
%         r = 6;
%     else
%         r = 1;
%     end
end

function dist = distance(refseq, seq)
% This function calculates the Hamming distance between the reference
% sequence and another sequence.
    dist = sum(refseq ~= seq);
end

function [seq, aseq, ntot, V, r, I, d] = ...
    replicate(i0, seq, aseq, mu, ntot, V, r, I, d, ksi, sigma, r0, refSeq)
% This function lets a virus of strain i0 replicate. During the
% replication, each nucleotide has a probability mu to mutate.
% If a new strain is formed by mutation, the fitness of this new strain is
% calculated with the fitness function.
    oldsequence = seq{i0};
    newsequence = oldsequence;
    L = length(oldsequence);
    nt = ['a','t','g','c'];
    % loop through nucleotides in sequence
    for n = 1:L
        % determine if a mutation occurs at this location:
        if rand < mu
            newnt = nt(nt ~= newsequence(n));
            newsequence(n) = newnt(randi(3));
        end
    end
    if strcmp(oldsequence, newsequence)
        V(i0) = V(i0) + 1;
        return
    end
    % determine if the mutated strain was already present
    sameSeq = find(strcmp(seq, newsequence));
    if ~isempty(sameSeq) 
        % if yes, then add a new free virus to that strain
        i0 = sameSeq;
        V(i0) = V(i0) + 1;
    else
        % if not, then create a new viral sequence and determine its
        % fitness
        ntot = ntot + 1;
        seq{ntot} = newsequence;
        aseq{ntot} = nt2aa(newsequence);
        V(ntot) = 1;
        di = distance(refSeq, aseq{ntot});
        d(ntot) = di;
        r(ntot) = fitness(di, r0, ksi, sigma);
        I(ntot) = 0;
    end
end

function [seq, aseq, V, I, r, d, ntot] = ...
             remove(seq, aseq, V, I, r, d, ntot)
    extinct = find(V+I == 0);  
    if ~isempty(extinct)
        V(extinct) = [];
        I(extinct) = [];
        seq(extinct) = [];
        aseq(extinct) =[];
        r(extinct) = [];
        d(extinct) = [];
        ntot = ntot - length(extinct);
    end
end
