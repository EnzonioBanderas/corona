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

clear all 

% Initialize ##############################################################
La = 3.0;       % length of AA sequence
a = 1.0e-3;     % infection rate per day.
b = 3.0;        % death/clearence rate of infected cells per day.
mu = 1.0e-4;    % mutation rate per nt per generation.
U = 1.0e5;      % initial number of uninfected cells.

V(1) = 400;     % initial number of free viral particles
I(1) = 0;       % Initial number of infected cells
r0 = 6.0;       % Fitness of WT sequence 
r(1) = r0;      % Array to keep track of fitnesses 

ntot = 1;       % initial number of distinct viral genotypes

ksi = 1.0;    	% fitness decay
sigma = 0.1;    % standard deviation of fitness

seq{1} = 'catacacaagga';    % initial genome sequence
aseq{1} = nt2aa(seq{1});    % initial aa sequence

T = 8;          % maximal time (days)
t = 0.0;        % initial time (days)

s = 1;          % counter
data(s, :) = [t, ntot, U, sum(I), sum(V)];

% Start loop ##############################################################
while t < T
     if sum(I+V) == 0                                                       % Stop if there are no more viral particles left
         break
     end
     s = s + 1;                                                             % counter
     % calculate total reaction rate
     Rtot = (a*U + b)*sum(V) + b*sum(I) + sum(r.*I);                        
     % take time step
     dt = (1/Rtot) * log(1/rand);                                           % time step is exponentially distributed and depends on total reaction rate
     t = t + dt;
     % choose a random number between 0 and Rtot
     rnd = Rtot * rand;
     c = 0;
     % loop over viral strains
     for i = 1:ntot    
         % choose which of the 4 reactions will occur and for which strain
         c = c + a*U*V(i);                                                  % R1: infection by strain i
         if c > rnd                                                         
             I(i) = I(i) + 1;                                               
             V(i) = V(i) - 1;
             U = U - 1;                                                    
             break
         end
         c = c + r(i)*I(i);                                                 % R2: replication of strain i
         if c > rnd
             i0 = i;
             [seq, aseq, mu, ntot, V, r, I] = ...
                 replicate(i0, seq, aseq, mu, ntot, V, r, I, ksi, sigma,r0);
             break
         end
         c = c + b*I(i);                                                    % R3: clearence of a cell infected by strain i
         if c > rnd
             I(i) = I(i) - 1;
             break
         end
         c = c + b*V(i);                                                    % R4: clearence of a free viral particle of strain i
         if c > rnd
             V(i) = V(i) - 1;
             break
         end
     end % for-loop
     % remove strains that have gone extinct
     extinct = find(V+I == 0);  
     if ~isempty(extinct)
         V(extinct) = [];
         I(extinct) = [];
         seq(extinct) = [];
         r(extinct) = [];
         ntot = ntot - length(extinct);
     end
     % collect statistics
     data(s,:) = [t, ntot, U, sum(I), sum(V)];
     % display progression
     if mod(s,1000) == 0
         disp(['t=',num2str(t),', ntot=',num2str(ntot)])
     end
end % while-loop
% plot the results ########################################################
figure()
plot(data(:,1), data(:,2))
xlabel('Time (days)')
ylabel('Number of viral strains')
title('Number of viral strains vs time')
figure()
semilogy(data(:,1), data(:,3), data(:,1), data(:,4))
xlabel('time (days)')
ylabel('Cell count')
legend('Number of uninfected cells', 'Number of infected cells')
title('Viral load vs time')

%% Functions 
function r = fitness(d, r0, ksi, sigma)
% This function determines the fitness of a strain with a distance d from
% the reference sequence, which has fitness r0.
% The fitness is normal distributed with mean mu = r0*exp(-d/ksi) and
% standard deviation sigma.
    mu = r0*exp(-d/ksi);                                        
    r = sigma*randn + mu;
end

function dist = distance(refseq, seq)
% This function calculates the Hamming distance between the reference
% sequence and another sequence.
    dist = sum(refseq ~= seq);
end

function [seq, aseq, mu, ntot, V, r, I] = ...
    replicate(i0, seq, aseq, mu, ntot, V, r, I, ksi, sigma, r0)
% This function lets a virus of strain i0 replicate. During the
% replication, each nucleotide has a probability mu to mutate.
% If a new strain is formed by mutation, the fitness of this new strain is
% calculated with the fitness function.
    gsequence = seq{i0};
    L = length(gsequence);
    nt = ['a','t','g','c'];
    % loop through nucleotides in sequence
    for n = 1:L
        % determine if a mutation occurs at this location:
        if rand < mu
            newnt = nt(nt ~= gsequence(n));
            gsequence(n) = newnt(randi(3));
        end
    end
    % determine if the mutated strain was already present
    sameSeq = find(strcmp(seq, gsequence));
    if ~isempty(sameSeq) 
        % if yes, then add a new free virus to that strain
        i0 = sameSeq;
        V(i0) = V(i0) + 1;
    else
        % if not, then create a new viral sequence and determine its
        % fitness
        ntot = ntot + 1;
        seq{ntot} = gsequence;
        V(ntot) = 1;
        d = distance(seq{1}, gsequence);
        r(ntot) = fitness(d, r0, ksi, sigma);
        I(ntot) = 0;
    end
end