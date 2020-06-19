clear all
% close all

%% Viral evolution with Gillespie algorithm

% Source of model: 
% Woo & Reifman (2013). Quantitative Modeling of Virus Evolutionary
% Dynamics and Adaptation in Serial Passages Using Empirically Inferred
% Fitness Landscapes. Journal of Virology, V88 N2, p1039 - 1050.

% Possible events/reactions:
% Infection: U + Vn -> In     (R1) rate = a, reaction_rate = a*U*Vn 
% Replication: In -> In + Vm  (R2) rate = r(d), reaction_rate = r*I
% Death : In -> 0             (R3) rate = b
% Clearence: Vn -> 0          (R4) rate = b
%% Get protein and genome reference sequences #############################
[gRefSeq, pRefSeq, pNames, pInfo, proteinLocation, genomeLocation] = getRefSeq();
translateCodon = geneticcode();
[beta0, beta1] = logisticRegressionProteins();
% Store betas of logistic regression in pInfo
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
U0 = 1e4;                   % initial number of not-infected cells
ksi = 1.0;                 	% fitness decay
sigma = 0.1;                % standard deviation of fitness
T = inf;                    % maximal time (days)
mu_array = linspace(1e-6, 1e-3,5);                 % mutation rate(s). Can be a single float (one mutation rate) or an array.
mu_array = 1e-3;
N_mu = length(mu_array);                            % number of mutation rates to test
d0 = zeros(length(pNames),1);                       % distance of WTseq to fittest strain
r0 = 2;                                             % Fitness of reference sequence 

% Error threshold detectors
Uarray = zeros(1, N_mu);                            % store the number of not-infected cells at the end of the infection
statD = zeros(1, N_mu);                             % stationary value of number of mutations
statR = zeros(1, N_mu);                             % stationary value of relative fitness
statY = zeros(N_mu, La + 1);                    	% stationary value of particle number
titer = struct();

S = 1e3;                                            % initial size of arrays

% Start for-loop over mutation rates ######################################
for k = 1:N_mu
    mu = mu_array(k);                                                       % mutation rate
    % initialize
    disp(mu)
    U = U0;                                                                 % number of infected cells.
    seq_loc = cell(1,S); aseq_loc = cell(1,S);                                % number of infected cells.
    seq_mut = cell(1,S); aseq_mut = cell(1,S);
    seq_nMut = zeros(1,S);
%     seq{1} = cell(1,2); aseq{1} = cell(1,2);
%     seq{1} = [];                                                     	% cell array to keep track of nt sequences
%     aseq{1} = [];                                                      % cell array to keep track of aa sequences   
    % counters
    ntot = 1;                                                               % initial number of distinct viral genotypes
    nAA = 1;                                                                % initial number of distinct viral phenotypes
    s = 0;                                                                  % counter for while loop
    m = 1;                                                                  % counter to collect stats
    t = 0.0;                                                                % initial time (days)
    % arrays to keep track of viral strains
    V = zeros(1, S);                                                        % number of free viral particles per strain
    V(1) = V0;
    I = zeros(1, S);                                                        % number of infected cells per strain
    r = zeros(1, S);                                                        % fitness per strain
    r(1) = r0;                                                              
    d = zeros(length(pNames), S);                                           % distance to reference per strain per protein
    d(:,1) = d0;
    dtot = zeros(1, S);                                                     % total distance to reference (sum of all proteins) per strain
    % arrays to collect stats
    Y = zeros(S, La + 1);                                                   % matrix for number of viruses with distance d from WT seq.
   	Y(1,1) = V(1);                                                      	
    
    data = zeros(S, 6);                                                     % matrix to collect time, # distinct genotypes, # distinct phenotypes, # not-infected cells, # infected cells, # viruses.
  	data(1,:) = [t, ntot, nAA, U, sum(I), sum(V)];                          
    meanDistance = zeros(1,S);                                              % mean number of mutations
    meanFitness = zeros(1,S);                                               % mean fitness of population
 	meanFitness(1) = r0; 
    maxD = zeros(1,S);
    maxR = zeros(1,S);
    maxR(1) = r0;
    
    % Start while loop ####################################################
    while t < T
         % Stop if there are no more viral particles left:
         if sum(I+V) == 0                                                  	
           	break
         end            
         s = s + 1;                                                         
         % calculate total reaction rate
         Rtot = (a*U + b)*sum(V) + b*sum(I) + sum(r.*I);                                          
         % take time step
         dt = (1/Rtot) * log(1/rand);                                       % time step is exponentially distributed and depends on total reaction rate
         t = t + dt;
         % choose a random number between 0 and Rtot
         rnd = Rtot * rand;
         c = 0;
         % Increase the size of arrays if it is necessary
         if find(V+I, 1,'last') + 10 > length(V)
             V = [V, zeros(1,S)];
             I = [I, zeros(1,S)];
             r = [r, zeros(1,S)];
             d = [d, zeros(length(pNames), S)];
             dtot = [dtot, zeros(1,S)];
             seq_loc = [seq_loc, cell(1,S)]; aseq_loc = [aseq_loc, cell(1,S)];
             seq_mut = [seq_mut, cell(1,S)]; aseq_mut = [aseq_mut, cell(1,S)];
             seq_nMut = [seq_nMut, zeros(1,S)];
         end
         % loop over viral strains
         for i = 1:find(V+I, 1,'last')    
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
                 [seq_loc, seq_mut, aseq_loc, aseq_mut, ntot, nAA, V, r, I, d, dtot, seq_nMut] = ...
                     replicate(i0, seq_loc, seq_mut, aseq_loc, aseq_mut, mu, ...
                     ntot, nAA, V, r, I, d, dtot, ... 
                     sigma,r0, gRefSeq, L, pRefSeq, pInfo, seq_nMut, proteinLocation, translateCodon);
                 break
             end
             c = c + b*I(i);                                                % R3: clearence of a cell infected by strain i
             if c > rnd
                 I(i) = I(i) - 1;
                 if V(i) + I(i) == 0    % if the strain has gone extinct, remove it
                     i0 = i;
                     [seq_loc, seq_mut, aseq_loc, aseq_mut, V, I, r, d, dtot, ntot, nAA, seq_nMut] = ...
                         remove(i0, seq_loc, seq_mut, aseq_loc, aseq_mut, V, I, r, d, dtot, ntot, nAA, seq_nMut);
                 end
                 break % for-loop
             end
             c = c + b*V(i);                                                % R4: clearence of a free viral particle of strain i
             if c > rnd
                 V(i) = V(i) - 1;
                 if V(i) + I(i) == 0    % if the strain has gone extinct, remove it
                     i0 = i;
                     [seq_loc, seq_mut, aseq_loc, aseq_mut, V, I, r, d, dtot, ntot, nAA, seq_nMut] = ...
                         remove(i0, seq_loc, seq_mut, aseq_loc, aseq_mut, V, I, r, d, dtot, ntot, nAA, seq_nMut);
                 end
                 break % for-loop
             end
         end % for-loop   
         % Collect statistics every 100 iterations:
         if mod(s,100) == 0
             m = m + 1;
             % Increase size of arrays if necessary
             if m > length(meanFitness)
                 Y = [Y; zeros(S, La + 1)];
                 data = [data; zeros(S, 6)];
                 meanDistance = [meanDistance, zeros(1,S)];
                 meanFitness = [meanFitness, zeros(1,S)];
                 maxD = [maxD, zeros(1,S)];
                 maxR = [maxR, zeros(1,S)];
             end
             % Group strains with the same distance and store them in Y:
             alive = (V + I ~= 0);
             shortV = V(alive); % remove zeros of extinct strains
             shortI = I(alive);
             for j = 1:max(dtot)+1
                dist = j - 1;
                d_logical = (dtot(alive) == dist);
                sumV = sum(shortV(d_logical));
                sumI = sum(shortI(d_logical));
                Y(m,j) = sumV + sumI;
             end
             % Number of distinct amino acid sequences:
             % nAA = length(unique(aseq(~cellfun('isempty', aseq))));  
             % Collect all data
             meanFitness(m) = sum(r.*(V+I)) / sum(V+I);
             meanDistance(m) = sum(dtot.*(V+I)) / sum(V+I);
             data(m,:) = [t, ntot, nAA, U, sum(I), sum(V)];
             maxD(m) = max(dtot);
             maxR(m) = max(r);
         end
         %display progression every 1000 iterations
         if mod(s,1000) == 0
             disp(['t=',num2str(t),', ntot=',num2str(ntot)])
         end
    end % while-loop
    % Remove zeros at the end of vectors
	Y = Y(1:m, :);
    data = data(1:m, :);
    meanFitness = meanFitness(1:m);
    meanDistance = meanDistance(1:m);
    maxD = maxD(1:m);
    maxR = maxR(1:m);
    % Calculate stationary values:
    time = data(:,1)';
    delta_t = time(2:end) - time(1:end-1);
    sumY = 1 ./ sum(Y,2);
    relativeY = Y .* repmat(sumY, 1, size(Y,2));             
    statY(k,:) = (delta_t * relativeY(2:end, :)) / time(end) ;
    statD(k) = sum(meanDistance(2:end) .* delta_t) / time(end);
    statR(k) = sum(meanFitness(2:end) .* delta_t) / time(end); 
    % Collect error threshold detectors
    Uarray(k) = U;
    titer(k).a = a;
    titer(k).Mu = mu;
   	titer(k).Time = data(:,1);
    titer(k).ntot = V;
    titer(k).I = I;
    titer(k).data = data;
    titer(k).Load = data(:,5) + data(:,6);
end % for-loop
save('Gillespie.mat','-v7.3');