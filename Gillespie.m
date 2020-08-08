% clear all
% close all

%function Gillespie(x)
x = '1';
myStream = RandStream('mlfg6331_64', 'Seed', str2num(x));

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
[gRefSeq, pRefSeq, pNames, proteinLocation, genomeLocation] = getRefSeq();
translateCodon = geneticcode();
[beta, sigma] = logisticRegressionProteins();

%% Initialize ############################################################# 
a = 2e-5;                   % infection rate per day.
b = 0.4;                    % death/clearance rate of infected cells per day.
V0 = 400;                   % initial number of viruses
U0 = 1e5;                   % initial number of uninfected cells
ksi = 1.0;                 	% fitness decay
T = inf;                    % maximal time (days)
t_anti = inf;
mu = 1e-6;
mu_array = mu;
%mu_array = linspace(1e-6, 1e-3,5);                 % mutation rate(s). Can be a single float (one mutation rate) or an array.
r0 = 2;                                             % Fitness of reference sequence 
alpha_array = 1;
antiviral = false;

S = 1e3;                    % initial size of arrays
nIter_record = 100;         % record every 100 iterations
nIter_print = 1000;         % print every 1000 iterations
kill = 0;                      % killing of infected cells by immune cells
stim = 0;                   % stimulation of immune cells by infected cells

% Follows
N_mu = length(mu_array);                            % number of mutation rates to test
L = length(gRefSeq);                                % length of genome sequence
La = length(pRefSeq);  
d0 = zeros(length(pNames),1);                       % distance of WTseq to fittest strain

% Error threshold detectors
Uarray = zeros(1, N_mu);                            % store the number of not-infected cells at the end of the infection
statD = zeros(1, N_mu);                             % stationary value of number of mutations
statR = zeros(1, N_mu);                             % stationary value of relative fitness
statY = zeros(N_mu, La + 1);                    	% stationary value of particle number

titer = struct();



% Start for-loop over mutation rates ######################################
for k = 1:N_mu
    r0 = 2
    %mu = mu_array(k);                                                       % mutation rate
    alpha = alpha_array(k)
    % initialize
    disp(mu)
    U = U0;                                                                 % number of infected cells.
    seq_loc = cell(1,S); aseq_loc = cell(1,S);                                % number of infected cells.
    seq_mut = cell(1,S); aseq_mut = cell(1,S);
    seq_nMut = zeros(1,S); aseq_nMut = zeros(1,S);                         % cell array to keep track of aa sequences   
    
    aseqUniq_loc = cell(1,S);
    aseqUniq_mut = cell(1,S);
    aseqUniq_nMut = zeros(1,S);
    aseqUniq_n = zeros(1,S); aseqUniq_n(1) = 1;
    aseqUniq_i = cell(1,S); aseqUniq_i{1} = 1;
    aseqUniq_r = zeros(1,S); aseqUniq_r(1) = r0;
    
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
    Tcells = zeros(1, S);
    Tcells_perStrain = zeros(1, S);
    r = zeros(1, S);                                                        % fitness per strain
    r(1) = r0;                                                              
    d = zeros(length(pNames), S);                                           % distance to reference per strain per protein
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
         if t >= t_anti && antiviral == false
             r0 = alpha*r0;
             r  = alpha*r;
             antiviral = true;
         end    
         s = s + 1;                                                         
         % calculate total reaction rate
         Rtot = (a*U + b)*sum(V) + b*sum(I) + sum(r.*I) + ...
                stim*sum(I) + b*sum(Tcells) + sum(k*Tcells_perStrain.*I); 
         % take time step
         dt = (1/Rtot) * log(1/rand(myStream));                                       % time step is exponentially distributed and depends on total reaction rate
         t = t + dt;
         % choose a random number between 0 and Rtot
         rnd = Rtot * rand(myStream);
         c = 0;
         % Increase the size of arrays if it is necessary
         existing = find(V+I);
         if existing(end) + 10 > length(V)
             V = [V, zeros(1,S)];
             I = [I, zeros(1,S)];
             Tcells_perStrain = [Tcells_perStrain, zeros(1, S)];
             r = [r, zeros(1,S)];
             d = [d, zeros(length(pNames), S)];
             dtot = [dtot, zeros(1,S)];
             seq_loc = [seq_loc, cell(1,S)]; aseq_loc = [aseq_loc, cell(1,S)];
             seq_mut = [seq_mut, cell(1,S)]; aseq_mut = [aseq_mut, cell(1,S)];
             seq_nMut = [seq_nMut, zeros(1,S)]; aseq_nMut = [aseq_nMut, zeros(1,S)];
         end
         existing_epitope = find(aseqUniq_n);
         if existing(end) + 10 > length(Tcells)
             Tcells = [Tcells, zeros(1,S)];
             aseqUniq_loc = [aseqUniq_loc, cell(1,S)];
             aseqUniq_mut = [aseqUniq_mut, cell(1,S)];
             aseqUniq_nMut = [aseqUniq_nMut, zeros(1,S)];
             aseqUniq_n = [aseqUniq_n, zeros(1,S)];
             aseqUniq_i = [aseqUniq_i, cell(1,S)];
             aseqUniq_r = [aseqUniq_r, zeros(1,S)];
         end
         % loop over viral strains
         reaction_happened = false;
         for i = existing    %1:find(V+I, 1,'last')?   
             % choose which of the 4 reactions will occur 
             % and for which strain
             c = c + a*U*V(i);                                              % R1: infection by strain i
             if c > rnd                                                         
                 I(i) = I(i) + 1;                                               
                 V(i) = V(i) - 1;
                 U = U - 1;      
                 reaction_happened = true;      
                 break
             end
             c = c + r(i)*I(i);                                             % R2: replication of strain i
             if c > rnd
                 i0 = i;
                [seq_loc, seq_mut, seq_nMut, ...
                    aseq_loc, aseq_mut, aseq_nMut, ...
                    ntot, nAA, V, r, I, d, dtot, ...
                    aseqUniq_loc, aseqUniq_mut, aseqUniq_nMut, ...
                    aseqUniq_n, aseqUniq_i, aseqUniq_r] = ...
                replicate(i0, myStream, ...
                    seq_loc, seq_mut, seq_nMut, ...
                    aseq_loc, aseq_mut, aseq_nMut, ...
                    mu, ntot, nAA, V, r, I, d, dtot, ...
                    sigma, r0, ...
                    gRefSeq, L, pRefSeq, beta, proteinLocation, translateCodon, ...
                    aseqUniq_loc, aseqUniq_mut, aseqUniq_nMut, ...
                    aseqUniq_n, aseqUniq_i, aseqUniq_r);
                 reaction_happened = true;
                break
             end
             c = c + b*I(i);                                                % R3: clearence of a cell infected by strain i
             if c > rnd
                 I(i) = I(i) - 1;
%                  disp(['V=',num2str(V(i)),', I=',num2str(I(i))])
                 if V(i) + I(i) == 0    % if the strain has gone extinct, remove it
                     i0 = i;
                    [seq_loc, seq_mut, seq_nMut, ...
                        aseq_loc, aseq_mut, aseq_nMut, ...
                        V, I, r, d, dtot, ntot, nAA, ...
                        Tcells, ...
                        aseqUniq_loc, aseqUniq_mut, aseqUniq_nMut, ...
                        aseqUniq_n, aseqUniq_i, aseqUniq_r] = ...
                    remove(i0, ...
                        seq_loc, seq_mut, seq_nMut, ...
                        aseq_loc, aseq_mut, aseq_nMut, ...
                        V, I, r, d, dtot, ntot, nAA, ...
                        Tcells, ...
                        aseqUniq_loc, aseqUniq_mut, aseqUniq_nMut, ...
                        aseqUniq_n, aseqUniq_i, aseqUniq_r);
                 end
                 reaction_happened = true;
                 break % for-loop
             end
             c = c + b*V(i);                                                % R4: clearence of a free viral particle of strain i
             if c > rnd
                 V(i) = V(i) - 1;
%                  disp(['V=',num2str(V(i)),', I=',num2str(I(i))])
                 if V(i) + I(i) == 0    % if the strain has gone extinct, remove it
                     i0 = i;
                    [seq_loc, seq_mut, seq_nMut, ...
                        aseq_loc, aseq_mut, aseq_nMut, ...
                        V, I, r, d, dtot, ntot, nAA, ...
                        Tcells, ...
                        aseqUniq_loc, aseqUniq_mut, aseqUniq_nMut, ...
                        aseqUniq_n, aseqUniq_i, aseqUniq_r] = ...
                    remove(i0, ...
                        seq_loc, seq_mut, seq_nMut, ...
                        aseq_loc, aseq_mut, aseq_nMut, ...
                        V, I, r, d, dtot, ntot, nAA, ...
                        Tcells, ...
                        aseqUniq_loc, aseqUniq_mut, aseqUniq_nMut, ...
                        aseqUniq_n, aseqUniq_i, aseqUniq_r);
                 end
                 reaction_happened = true;
                 break % for-loop
             end
         end % for-loop 
         if ~reaction_happened
         for iEpitope = existing_epitope % loop over unique epitopes 
             Epitope2Strain_index = aseqUniq_i{iEpitope};
             Epitope2Strain_index = Epitope2Strain_index(I(Epitope2Strain_index)~=0);
             c = c + stim*sum(I(Epitope2Strain_index));                                             % Stimulation: Infected cell stimulates formation of T-cell specific to virus
             if c > rnd
                 Tcells(iEpitope) = Tcells(iEpitope) + 1;
                 Tcells_perStrain(Epitope2Strain_index) = Tcells_perStrain(Epitope2Strain_index) + 1;
                 break
             end 
             c = c + kill*sum(I(Epitope2Strain_index))*Tcells(iEpitope);                                           % Killing: Cytotoxic T cells kill infected cells
             if c > rnd
                 i0 = Epitope2Strain_index(randi(myStream, length(Epitope2Strain_index)));
                 I(i0) = I(i0) - 1;
                 if V(i0) + I(i0) == 0    % if the strain has gone extinct, remove it
                    [seq_loc, seq_mut, seq_nMut, ...
                        aseq_loc, aseq_mut, aseq_nMut, ...
                        V, I, r, d, dtot, ntot, nAA, ...
                        Tcells, ...
                        aseqUniq_loc, aseqUniq_mut, aseqUniq_nMut, ...
                        aseqUniq_n, aseqUniq_i, aseqUniq_r] = ...
                    remove(i0, ...
                        seq_loc, seq_mut, seq_nMut, ...
                        aseq_loc, aseq_mut, aseq_nMut, ...
                        V, I, r, d, dtot, ntot, nAA, ...
                        Tcells, ...
                        aseqUniq_loc, aseqUniq_mut, aseqUniq_nMut, ...
                        aseqUniq_n, aseqUniq_i, aseqUniq_r);
                 end
                 break
             end
             c = c + b*Tcells(iEpitope);                                                % Clearance: Immune cells die
             if c > rnd
                 Tcells(iEpitope) = Tcells(iEpitope) - 1;
                 Tcells_perStrain(Epitope2Strain_index) = Tcells_perStrain(Epitope2Strain_index) - 1;
%                Use custom remove_immune instead to keep AAs in memory as
%                long as Tcells exist that are specific to those AAs
                 break
             end
         end
         end
         
         % Collect statistics every 100 iterations:
         if mod(s, nIter_record) == 0
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
         if mod(s, nIter_print) == 0
             disp(['t=',num2str(t),', \t ntot=',num2str(ntot), ', \t U=', num2str(U)])
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
    titer(k).alpha = alpha;
    titer(k).Mu = mu;
    titer(k).t = t;
    titer(k).ntot = ntot;
    titer(k).nAA = nAA;
    titer(k).U_sum = sum(U);
    titer(k).I_sum = sum(I);
    titer(k).V_sum = sum(V);
    titer(k).data = data;
    titer(k).statY = statY(k);
    titer(k).statD = statD(k);
    titer(k).statR = statR(k);
    titer(k).maxD = maxD;
    titer(k).maxR = maxR;
end % for-loop
% fname = ['Gillespie', num2str(x), '.mat']
% save(fname, 'titer','mu','statY','statD','statR','maxD','maxR',...      
%     '-v7.3');
fname = ['Output/Gillespie.mat']
save(fname, 'titer','mu','statY','statD','statR','maxD','maxR',...      
    '-v7.3');
%end 