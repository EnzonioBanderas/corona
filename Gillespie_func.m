function Gillespie_func(x)
    
%     x = '1';
    disp(x)
%     nPoint_prev = 1e6;
    myStream = RandStream('mlfg6331_64', 'Seed', str2num(x));
%     myStream = RandStream('mlfg6331_64', 'Seed', str2num(x)+nPoint_prev);

%     p = sobolset(1, 'Skip', str2num(x)-1);
%     p = scramble(p,'MatousekAffineOwen');
%     rand_sobol = net(p, 1); % generate 1 points of the sobol sequence
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Data', filesep, 'params_lhs_1e4_7_Iter1e3.mat'], 'params_lhs') % load params_lhs
%     params_lhs = lhsdesign(1e3, 1, 'Iterations', 1e4); % [V0 or mu]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nLHS = size(params_lhs, 1);
    iLHS = mod(str2num(x)-1,nLHS)+1;

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
    %% Get protein and genome reference sequences 
    [gRefSeq, pRefSeq, pNames, proteinLocation, genomeLocation] = getRefSeq();
    translateCodon = geneticcode();
    [beta, sigma] = logisticRegressionProteins();

    %% Initialize 
%     U0 = 1e6;                   % initial number of uninfected cells
    distribution = 'normal';
    
    b_center = 0.9;
    r0_center = 1.5;    
    mu_center = 1e-6;
    U0_center = 1e4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if U0_center == 1e4
        a_center =  4.5e-3 ;
        V0_center = 40;
    elseif U0_center == 1e5
        a_center =  4.5e-4 ;
        V0_center = 400;
    elseif U0_center == 1e6
        a_center =  4.5e-5 ;
        V0_center = 4e3;
    end
    alpha = 1;
    
%     a = a_center;
%     a = 10^((rand_sobol*-3)-2); % 10^-2 to 10^-4 (log)
%     a = (rand_sobol*6e-3)+1.5e-3; % 3*10^-3 to 6*10^-3 (lin)
%     a = (params_lhs(iLHS, 1)*8e-3)+1e-3; % 1*10^-3 to 9*10^-3 (lhs lin)
    a = (params_lhs(iLHS, 1)*a_center*0.4)+a_center*0.8; % 1*10^-3 to 9*10^-3 (lhs lin loc)
%     a = (params_lhs(iLHS, 1)*a_center*2.9)+a_center*0.1; % 1*10^-3 to 9*10^-3 (lhs lin glob)
    
%     b = b_center;         
%     b = 10^((rand_sobol*3)-2); % 1e-2 to 1e1 (log)
%     b = (rand_sobol*2.9)+0.1; % 0.1 to 3 (lin)
%     b = (params_lhs(iLHS, 2)*2.9)+0.1; % 0.1 to 3 (lhs lin)
    b = (params_lhs(iLHS, 2)*b_center*0.4)+b_center*0.8; % 1*10^-3 to 9*10^-3 (lhs lin loc)
%     b = (params_lhs(iLHS, 2)*b_center*2.9)+b_center*0.1; % 1*10^-3 to 9*10^-3 (lhs lin glob)

%     c = b;  
    c_center = b_center;
    c = (params_lhs(iLHS, 7)*c_center*0.4)+c_center*0.8; % 8e3 to 12e3 (lhs lin loc)
    
%     r0 = r0_center;
%     r0 = 10^((rand_sobol*2)-1); % 1e-1 to 1e1 (log)
%     r0 = (rand_sobol*2.9)+0.1; % 0.1 to 3 (lin)
%     r0 = (params_lhs(iLHS, 3)*2.9)+1; % 1 to 3 (lhs lin) 
    r0 = (params_lhs(iLHS, 3)*r0_center*0.4)+r0_center*0.8; % 1*10^-3 to 9*10^-3 (lhs lin loc)
%     r0 = (params_lhs(iLHS, 3)*r0_center*2.9)+r0_center*0.1; % 1*10^-3 to 9*10^-3 (lhs lin glob)

%     mu = mu_center;
%     mu = 10^((params_lhs(iLHS)*-6)-2); % 1e-8 to 1e-2 (log)
%     mu = (rand_sobol*1e-6)+0.5e-6; % 0.5e-6 to 1.5e-6 (lin)
%     mu = (params_lhs(iLHS, 4)*2.9e-6)+0.1e-6; % 0.1e-6 to 3.0e-6 (lhs lin)
    mu = (params_lhs(iLHS, 4)*mu_center*0.4)+mu_center*0.8; % 1*10^-3 to 9*10^-3 (lhs lin loc)
%     mu = (params_lhs(iLHS, 4)*mu_center*2.9)+mu_center*0.1; % 1*10^-3 to 9*10^-3 (lhs lin glob)

%     V0 = V0_center;
%     V0 = round(4*(10^((rand_sobol*4)+0))); % 4e0 to 4e4 (log)
%     V0 = round((rand_sobol*90)+10); % 10 to 100 (lin)
    V0 = round((params_lhs(iLHS, 5)*V0_center*0.4)+V0_center*0.8); % 32 to 48 (lhs lin loc)
    
%     U0 = U0_center;
%     U0 = 10^((params_lhs(iLHS)*5)+2); % 1e2 to 1e7 (log)
    U0 = round((params_lhs(iLHS, 6)*U0_center*0.4)+U0_center*0.8); % 8e3 to 12e3 (lhs lin loc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    params = struct();
    params.('U0') = U0;
    params.('r0') = r0;
    params.('a') = a;
    params.('b') = b;
    params.('c') = c;
    params.('V0') = V0;
    params.('mu') = mu;

    T = 365*200;                    % maximal time (days) 200 years
    t_anti = inf;

    antiviral = false;

    kill = 0;                   % killing of infected cells by immune cells
    stim = 0;                   % stimulation of immune cells by infected cells

    % Follows
    nK = 1;                                            % Number of iterations to do with the same parameter values                    
    L = length(gRefSeq);                                % length of genome sequence
    La = length(pRefSeq);  
    d0 = zeros(length(pNames),1);                       % distance of WTseq to fittest strain
    data = struct();
    S = 1e3;                                            % initial size of arrays

    
    
    % Start for-loop over mutation rates 
    for k = 1:nK
        tic
        
        fprintf('\n Iteration %d: U0=%d \t V0=%d \t a=%f \t b=%f \t c=%f \t r0=%f \t mu=%e \n', k, U0, V0, a, b, c, r0, mu)
        
        % initialize
        U = U0;                                                                % number of infected cells.
        seq_loc = cell(1,S); aseq_loc = cell(1,S);                             % number of infected cells.
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
        t_collect = 0.0;
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

        data_collect = zeros(S, 11);                                                     % matrix to collect time, # distinct genotypes, # distinct phenotypes, # not-infected cells, # infected cells, # viruses. # infection occurrences # cell death occurences # virus clearance occurences # replication occurrences
        meanDistance = zeros(1,S);                                              % mean number of mutations
        meanFitness = zeros(1,S);                                               % mean fitness of population
        meanFitness(1) = r0; 
        maxD = zeros(1,S);
        maxR = zeros(1,S);
        maxR(1) = r0;
        ntot_e = zeros(1,S);
        ntot_e_I = zeros(1,S);
        aN = 0;
        bN = 0;
        cN = 0;
        rN = 0;
        data_collect(1,:) = [t, ntot, nAA, U, sum(I), sum(V), aN, bN, cN, rN, 0]; 
        diversity = zeros(1,S);

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
                    stim*sum(I) + b*sum(Tcells) + sum(kill*Tcells_perStrain.*I); 
             % take time step
             dt = (1/Rtot) * log(1/rand(myStream));                                       % time step is exponentially distributed and depends on total reaction rate
             t = t + dt;
             % choose a random number between 0 and Rtot
             rnd = Rtot * rand(myStream);
             z = 0;
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
                 z = z + a*U*V(i);                                              % R1: infection by strain i
                 if z > rnd                                                         
                     I(i) = I(i) + 1;                                               
                     V(i) = V(i) - 1;
                     U = U - 1;      
                     if U<0
                         fprintf('U=%d \t z=%e \t rnd=%e', U,z,rnd)
                     end
                     reaction_happened = true;      
                     aN = aN + 1;
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
                    rN = rN + 1;
                     reaction_happened = true;
                    break
                 end
                 z = z + c*I(i);                                                % R3: clearance of a cell infected by strain i
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
                     cN = cN + 1;
                     reaction_happened = true;
                     break % for-loop
                 end
                 z = z + b*V(i);                                                % R4: clearance of a free viral particle of strain i
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
                     bN = bN + 1;
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
                 if m > length(meanFitness)
                     Y = [Y; zeros(S, La + 1)];
                     data_collect = [data_collect; zeros(S, 11)];
                     meanDistance = [meanDistance, zeros(1,S)];
                     meanFitness = [meanFitness, zeros(1,S)];
                     maxD = [maxD, zeros(1,S)];
                     maxR = [maxR, zeros(1,S)];
                     ntot_e = [ntot_e, zeros(1,S)];
                     ntot_e_I = [ntot_e_I, zeros(1,S)];
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

                 meanFitness(m) = sum(r(alive) .* viralTiter) / sum(viralTiter);
                 meanDistance(m) = sum(dtot(alive) .* viralTiter) / sum(viralTiter);
                 data_collect(m,:) = [t, ntot, nAA, U, sum(I), sum(V), aN, bN, cN, rN, sum(viralTiter)];
                 maxD(m) = max(dtot);
                 maxR(m) = max(r);

                 diversity(m) = 1 - sum((viralTiter/sum(viralTiter)).^2); % Simpson's index as diversity
                 
                 % V evenness by NT
                 ntot_e(m) = 1 - sum((V/sum(V)).^2); % Simpson's index as evenness
                 % I evenness by NT
                 ntot_e_I(m) = 1 - sum((I/sum(I)).^2); % Simpson's index as evenness
                 % V evenness by AA
                 % I evenness by AA

                 fprintf(['t=',num2str(t),', \t ntot=',num2str(ntot), ', \t U=', num2str(U), ' \n'])
             end

        end % while-loop

        % Remove zeros at the end of vectors
        maxY = find(sum(Y) > 0, 1, 'last');

        Y = Y(1:m, :);         % discard the last row of zeros   
        data_collect = data_collect(1:m, :);
        meanFitness = meanFitness(1:m);
        meanDistance = meanDistance(1:m);
        maxD = maxD(1:m);
        maxR = maxR(1:m);
        ntot_e = ntot_e(1:m);
        ntot_e_I = ntot_e_I(1:m);
        diversity = diversity(1:m);

        % Calculate stationary values:
        time = data_collect(:, 1)';
        delta_t = time(2:end) - time(1:end-1);
        sumY = 1 ./ sum(Y,2);
        relativeY = Y .* repmat(sumY, 1, size(Y,2));       
%         statY = sum(delta_t .* relativeY(2:end,:), 1) / time(end);
        statD = sum(meanDistance(2:end) .* delta_t) / time(end);
        statR = sum(meanFitness(2:end) .* delta_t) / time(end); 
        statDiv = sum(diversity(2:end) .* delta_t) / time(end);
        
        
        
        % Collect information for each simulation
        data(k).alpha = alpha;

        data(k).t_end = t;
        data(k).ntot_end = ntot;
        data(k).nAA_end = nAA;
        data(k).U_end = U;
        data(k).I_sum_end = sum(I);
        data(k).V_sum_end = sum(V);

        data(k).t = data_collect(:,1);
        data(k).ntot = data_collect(:,2);
        data(k).nAA = data_collect(:,3);
        data(k).U = data_collect(:,4);
        data(k).I_sum = data_collect(:,5);
        data(k).V_sum = data_collect(:,6);
        data(k).aN = data_collect(:,7);
        data(k).bN = data_collect(:,8);
        data(k).cN = data_collect(:,9);
        data(k).rN = data_collect(:,10);
        data(k).viralTiter = data_collect(:,11);

%         data(k).statY = statY;
        data(k).relativeY = relativeY(:, 1:maxY);
        data(k).statD = statD;
        data(k).statR = statR;
        data(k).statDiv = statDiv;
        data(k).maxD = maxD;
        data(k).maxR = maxR;
        data(k).diversity = diversity;
        data(k).evenness = ntot_e;
        data(k).evenness_I = ntot_e_I;
        [data(k).V_peak, peakIndex] = max(data(k).V_sum);
        data(k).V_peakTime = data(k).t(peakIndex);
        data(k).toc = toc;
    end % for-loop

    disp(['Sum of simulation times (s): ', num2str(sum([data.toc]))])
    
    simString = 'iter1';
%     fname = 'Output/Gillespie.mat';
%     save(fname, 'data', 'params', '-v7.3');
    mkdir(simString)
    save([simString, filesep, simString, '_', x, '.mat'], 'data', 'params', '-v7.3');
end 



