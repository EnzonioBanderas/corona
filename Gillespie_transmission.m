function Gillespie_transmission(x)

tic;
myStream = RandStream('mlfg6331_64', 'Seed', str2num(x));
outputFolder = ['Output_', x];
mkdir(outputFolder);

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
U0 = 1e5;                   % initial number of uninfected cells
distribution = 'normal';
Nind = 100;

if U0 == 1e5
    r0 = 1.5 ;
    a =  4.5e-04 ;
    b =  0.9 ;
    c =  0.9 ;
    V0 = 400;
elseif U0 == 1e4
    r0 = 1.5 ;
    a =  4.5e-3 ;
    b =  0.9 ;
    c =  0.9 ;
    V0 = 40;
 elseif U0 == 1e6
    r0 = 1.5 ;
    a =  4.5e-5 ;
    b =  0.9 ;
    c =  0.9 ;
    V0 = 4e3;
    
end
alpha = 1;
mu = 1e-6;

params = struct();
params.('U0') = U0;
params.('r0') = r0;
params.('a') = a;
params.('b') = b;
params.('c') = c;
params.('V0') = V0;
params.('mu') = mu;

T = inf;                    % maximal time (days)
t_anti = inf;
antiviral = false;

kill = 0;                   % killing of infected cells by immune cells
stim = 0;                   % stimulation of immune cells by infected cells

% Follows
L = length(gRefSeq);                                % length of genome sequence
La = length(pRefSeq);  
d0 = zeros(length(pNames),1);                       % distance of WTseq to fittest strain
S = 1e3;

titer = struct();

t_L = 4.6;                                          % Latent period (days)
t_I = 5;                                            % infectious period (days)

% arrays to keep track of viral strains
transV = zeros(1, S);                                                        % number of free viral particles per strain
transV(1) = V0;
transR = zeros(1, S);                                                        % fitness per strain
transR(1) = r0;                                                              
transD = zeros(length(pNames), S);                                           % distance to reference per strain per protein
transD(:,1) = d0;
transDtot = zeros(1, S);                                                     % total distance to reference (sum of all proteins) per strain
% arrays to collect stats
transY = zeros(S, La + 1);                                                   % matrix for number of viruses with distance d from WT seq.
transY(1,1) = V0;

transTcells = zeros(1, S);
transTcells_perStrain = zeros(1, S);

transSeq_loc = cell(1,S);
transAseq_loc = cell(1,S);
transSeq_mut = cell(1,S);
transAseq_mut = cell(1,S);
transSeq_nMut = zeros(1,S);        
transAseq_nMut = zeros(1,S);

transAseqUniq_loc = cell(1,S);
transAseqUniq_mut = cell(1,S);
transAseqUniq_nMut = zeros(1,S);
transAseqUniq_n = zeros(1,S); transAseqUniq_n(1) = 1;
transAseqUniq_i = cell(1,S); transAseqUniq_i{1} = 1;
transAseqUniq_r = zeros(1,S); transAseqUniq_r(1) = r0;

% Start for-loop over mutation rates ######################################
for k = 1:Nind
    % initialize
    disp(k)
    
    U = U0;                                                                 % number of infected cells.
    t = 0.0;                                                                % initial time (days)
    t_collect = t;

    % Arrays to keep track of viral strains
    I = zeros(1, S);                                                        % number of infected cells per strain
    V = transV;
    Y = transY;
    r = transR;
    d = transD;
    dtot = transDtot;

    seq_loc = transSeq_loc;
    aseq_loc = transAseq_loc;
    seq_mut = transSeq_mut;
    aseq_mut = transAseq_mut;
    seq_nMut = transSeq_nMut;
    aseq_nMut = transAseq_nMut;
    
    aseqUniq_loc = transAseqUniq_loc;
    aseqUniq_mut = transAseqUniq_mut;
    aseqUniq_nMut = transAseqUniq_nMut;
    aseqUniq_n = transAseqUniq_n;
    aseqUniq_i = transAseqUniq_i;
    aseqUniq_r = transAseqUniq_r;

    Tcells = transTcells;
    Tcells_perStrain = transTcells_perStrain;
    
    % counters
    ntot = sum(V + I > 0);                                              	% initial number of distinct viral genotypes
    nAA = sum(aseqUniq_n > 0);                                              % initial number of distinct viral phenotypes
    s = 0;                                                                  % counter for while loop
    m = 1;                                                                  % counter to collect stats

    % arrays to collect stats
    data = zeros(S, 6);                                                     % matrix to collect time, # distinct genotypes, # distinct phenotypes, # not-infected cells, # infected cells, # viruses.
    data(1,:) = [t, ntot, nAA, U, sum(I), sum(V)];     
    meanD = zeros(1,S);                                              % mean number of mutations
    meanD(1) = sum(dtot .* V) / sum(V);
    meanR = zeros(1,S);                                               % mean fitness of population
    meanR(1) = sum(r.* V) / sum(V);
    maxD = zeros(1,S);
    maxD(1) = max(dtot);
    maxR = zeros(1,S);
    maxR(1) = max(r);
    diversity = 1 - sum((V/sum(V)).^2);

    t_transmission = t_L + t_I * rand;                                     % transmission occurs at a random time between t_L and t_L+t_I
    transmissionBoolean = false; 

    endInfection = struct();
    transmission = struct();
        
    % Start while loop ####################################################
    while t < T
         if t >= t_anti && antiviral == false
             r0 = alpha*r0;
             r  = alpha*r;
             antiviral = true;
         end    
         s = s + 1;                                                         
         % calculate total reaction rate
         Rtot = (a*U + b)*sum(V) + c*sum(I) + sum(r.*I) + ...
                stim*sum(I) + b*sum(Tcells) + sum(k*Tcells_perStrain.*I); 
         % take time step
         dt = (1/Rtot) * log(1/rand(myStream));                                       % time step is exponentially distributed and depends on total reaction rate
         t = t + dt;
         % choose a random number between 0 and Rtot
         rnd = Rtot * rand(myStream);
         z = 0;
         % Increase the size of arrays if it is necessary
         existing = find(V+I);
         existing_epitope = find(aseqUniq_n);
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
             z = z + a*U*V(i);                                              % R1: infection by strain i
             if z > rnd                                                         
                 I(i) = I(i) + 1;                                               
                 V(i) = V(i) - 1;
                 U = U - 1;      
                 reaction_happened = true;      
                 break
             end
             z = z + r(i)*I(i);                                             % R2: replication of strain i
             if z > rnd
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
                    aseqUniq_n, aseqUniq_i, aseqUniq_r, distribution);
                 reaction_happened = true;
                break
             end
             z = z + c*I(i);                                                % R3: clearence of a cell infected by strain i
             if z > rnd
                 I(i) = I(i) - 1;
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
             z = z + b*V(i);                                                % R4: clearence of a free viral particle of strain i
             if z > rnd
                 V(i) = V(i) - 1;
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
             
             z = z + stim*sum(I(Epitope2Strain_index));                                             % Stimulation: Infected cell stimulates formation of T-cell specific to virus
             if z > rnd
                 Tcells(iEpitope) = Tcells(iEpitope) + 1;
                 Tcells_perStrain(Epitope2Strain_index) = Tcells_perStrain(Epitope2Strain_index) + 1;
                 break
             end 
             
             z = z + kill*sum(I(Epitope2Strain_index))*Tcells(iEpitope);                                           % Killing: Cytotoxic T cells kill infected cells
             if z > rnd
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
             
             z = z + b*Tcells(iEpitope);                                                % Clearance: Immune cells die
             if z > rnd
                 Tcells(iEpitope) = Tcells(iEpitope) - 1;
                 Tcells_perStrain(Epitope2Strain_index) = Tcells_perStrain(Epitope2Strain_index) - 1;
%                Use custom remove_immune instead to keep AAs in memory as
%                long as Tcells exist that are specific to those AAs
                 break
             end
         end
         end
         
        % Stop if there are no more viral particles left:
         if sum(I+V) == 0                                                  	
           	break
         end
         
         % Collect statistics every 0.1 days:
         if t - t_collect > 0.1
             t_collect = t;
             m = m + 1;
             % Increase size of arrays if necessary
             if m > length(meanR)
                 Y = [Y; zeros(S, La + 1)];
                 data = [data; zeros(S, 6)];
                 meanD = [meanD, zeros(1,S)];
                 meanR = [meanR, zeros(1,S)];
                 maxD = [maxD, zeros(1,S)];
                 maxR = [maxR, zeros(1,S)];
                 diversity = [diversity, zeros(1,S)];
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
     
             % Collect all data
             viralTiter = shortV + shortI;
             
             meanR(m) = sum(r(alive) .* viralTiter) / sum(viralTiter);
             meanD(m) = sum(dtot(alive) .* viralTiter) / sum(viralTiter);
             data(m,:) = [t, ntot, nAA, U, sum(I), sum(viralTiter)];
             maxD(m) = max(dtot);
             maxR(m) = max(r);

             diversity(m) = 1 - sum((viralTiter/sum(viralTiter)).^2); % Simpson's index as diversity
             
             %disp(['t=',num2str(t),', \t ntot=',num2str(ntot), ', \t U=', num2str(U)])
         end
         
          % transmission:
         if t >= t_transmission && transmissionBoolean == false
            transmissionBoolean = true;
            
            % get last occupied entry in V and I vectors 
            
            tempY = Y(1:m, :);

            % sample a random pool of 400 viruses.
            totV = sum(V);
            pool = randi(totV, [1, V0]);
            % initiate arrays for transmission
            transV = zeros(1, S); 
            transR = zeros(1, S);
            transD = zeros(length(pNames), S);
            transDtot = zeros(1, S);

            transSeq_loc = cell(1,S);
            transAseq_loc = cell(1,S);
            transSeq_mut = cell(1,S);
            transAseq_mut = cell(1,S);
            transSeq_nMut = zeros(1,S);
            transAseq_nMut = zeros(1,S);
            
            transAseqUniq_loc = cell(1,S);
            % add 'dummy' 1 values:
            transAseqUniq_mut = cell(1,S); transAseqUniq_mut(:) = {'1'};
            transAseqUniq_nMut = zeros(1,S);
            transAseqUniq_n = zeros(1,S); 
            transAseqUniq_i = cell(1,S); 
            transAseqUniq_r = zeros(1,S);

            g = 0;
            p = 0;
            nu = 1;
            for i = 1:length(V)
                newg = g + V(i);
                numTransmitted = sum(pool > g & pool <= newg);
                if numTransmitted > 0
                    p = p + 1;
                    transV(p) = numTransmitted;
                 	transR(p) = r(i);
                    transD(:,p) = d(:,i);
                    transDtot(p) = dtot(i);

                    transSeq_loc{p} = seq_loc{i};
                    transSeq_mut{p} = seq_mut{i};
                    transSeq_nMut(p) = seq_nMut(i);
                    
                    transAseq_nMut(p) = aseq_nMut(i);
                    transAseq_loc{p} = aseq_loc{i};
                  	transAseq_mut{p} = aseq_mut{i};
                    
                    % Unique arrays:
                    alreadyPresent = false;
                    for j = 1 : nu
                        if isequal(transAseqUniq_loc{j}, aseq_loc{i})
                            if isequal(transAseqUniq_mut{j}, aseq_mut{i})
                            	alreadyPresent = true;
                                transAseqUniq_n(j) = transAseqUniq_n(j) + numTransmitted;
                                transAseqUniq_i{j} = [transAseqUniq_i{j}, p];
                                break % for loop over unique AA seqs
                            end  
                        end
                    end
                        
                    if ~alreadyPresent
                        transAseqUniq_loc{nu} = aseq_loc{i};
                        transAseqUniq_mut{nu} = aseq_mut{i};
                        transAseqUniq_nMut(nu) = aseq_nMut(i);
                        transAseqUniq_n(nu) = numTransmitted;
                        transAseqUniq_i{nu} = p;
                        transAseqUniq_r(nu) = r(i);
                        nu = nu + 1;
                    end                    
                end
                g = newg;
            end
            % remove the 'dummy' 1 values again
            transAseqUniq_mut(strcmp(transAseqUniq_mut, '1')) = {[]};
            
            % transmission Y array:
            alive = (transV ~= 0);
            shortV = transV(alive);
            transY = zeros(S, La + 1);                                     	% number of free viral particles per strain
            for j = 1:max(transDtot) + 1
                dist = j - 1;
                d_logical = (transDtot(alive) == dist);
                sumV = sum(shortV(d_logical));
                transY(1,j) = sumV;
            end      
            maxV = find(transV, 1,'last');

            transmission.('V') = transV(1:maxV);
            transmission.('r') = transR(1:maxV);
            transmission.('seq_loc') = transSeq_loc(1:maxV);
            transmission.('seq_mut') = transSeq_mut(1:maxV);
            transmission.('tTransmission') = t_transmission;
            transmission.('transD') = transD(:, 1:maxV);
            transmission.('transR') = transR(1:maxV);
            
            %transmission.('aseq_loc') = transAseq_loc(1:maxV);
            %transmission.('aseq_mut') = transAseq_mut(1:maxV);
        end % if transmission
             
    end % while-loop
    % Remove zeros at the end of vectors
    maxY = find(sum(Y) > 0, 1, 'last');
        
	Y = Y(1:m, :);         % discard the last row of zeros   
    data = data(1:m, :);
    meanR = meanR(1:m);
    meanD = meanD(1:m);
    maxD = maxD(1:m);
    maxR = maxR(1:m);
    diversity = diversity(1:m);
    
    % Calculate stationary values:
    time = data(1:m, 1)';
    delta_t = time(2:end) - time(1:end-1);
    sumY = 1 ./ sum(Y,2);
    relativeY = Y .* repmat(sumY, 1, size(Y,2));             
    statY = (delta_t * relativeY(2:end,:)) / time(end);
    statD = sum(meanD(2:end) .* delta_t) / time(end);
    statR = sum(meanR(2:end) .* delta_t) / time(end); 
    statDiv = sum(diversity(2:end) .* delta_t) / time(end);
    
    % Collect information for each simulation
    titer.alpha = alpha;
    titer.Mu = mu;
    titer.t = t;
    titer.ntot = ntot;
    titer.nAA = nAA;
    titer.U_sum = sum(U);
    titer.I_sum = sum(I);
    titer.V_sum = sum(V);
    titer.data = data;
    titer.statY = statY(1:maxY);
    titer.relativeY = relativeY(:, 1:maxY);
    titer.statD = statD;
    titer.statR = statR;
    titer.maxD = maxD;
    titer.maxR = maxR;
    titer.diversity = diversity;
    titer.statDiv = statDiv;
    
    fname = [outputFolder ,'/Gillespie', num2str(k), '.mat'];
    save(fname, 'titer', 'params', 'transmission', '-v7.3');
end % for-loop

disp(['Simulation time ', num2str(toc)])

end


