%{
Analysis of sequence and protein alignments.
%}

%% Load data
clear all
close all
clc

file = 'Sequence_alignment.mat';
data = load(file);

seqs = data.seqs;
names = data.names;

%% Find mutations
for s = 1:length(seqs)
    alignment = seqs(s).Alignment;
    Nmut = 0;
    bases = [];
    locations = [];
    for i = 1:size(alignment,2)                                             % loop through sequence alignment
        if alignment(1,i) ~= alignment(3,i) & alignment(1,i) ~= '-' ...     % check if there is a POINT mutation 
                & alignment(3,i) ~= 'N' & alignment(3,i) ~= '-' ...         % (indels and non-resolved nuclotides are not taken into account)
                & alignment(3,i) ~= 'R' & alignment(3,i) ~= 'Y' ...
                & alignment(3,i) ~= 'W' & alignment(3,i) ~= 'S' ...
                & alignment(3,i) ~= 'M' & alignment(3,i) ~= 'K' ...
                & alignment(3,i) ~= 'H' & alignment(3,i) ~= 'B' ...
                & alignment(3,i) ~= 'V' & alignment(3,i) ~= 'D'
            Nmut = Nmut+1;                                                  % number of mutations
            locations(Nmut) = i;                                            % location of mutations
            bases(:,Nmut) = alignment(:,i);                                 % base substitutions
        end
    end
    mutations(s) = Nmut;                                                    % list with all mutations
    bases = char(bases);                                                    
    seqs(s).('Mutations') = ...                                             % store information about mutations in the seqs structure.
        struct('Nm',Nmut,'Bases',bases, 'Locations',locations);
end
%% Display mutations in a histogram, and fit a distribution
figure()
histfit(mutations,max(mutations),'tlocations')                              % histfit, Nbins = max(mutations)
pd = fitdist(mutations', 'tlocations');
[h, p] = chi2gof(mutations','CDF',pd);
xlabel('Number of mutations')   
ylabel('Count')
title('Histogram of mutation number')

%% Store the mutation locations in an array
maxNm = max(mutations);                                                     % Maximal number of mutations
mutationArray = NaN(length(seqs), maxNm);                                   % rows in the array are viruses
for i=1:length(seqs)
    mutationVector = seqs(i).Mutations.Locations;                           % vector with mutation locations vor virus i
    mutationVector = [mutationVector, ...                                   % make the vector the right length
        NaN(1,maxNm-length(mutationVector))];                             
    mutationArray(i,:) = mutationVector;                                    % put vector in the mutation array
end
%% Plot the frequent mutations
thres = 100;                                                                % Enter a threshold
[frequentMutations, frequentMutationsCount] = ...                           % Find mutations that occur more frequent than threshold.     
    evaluateFrequentMutations(mutationArray, thres);                                                    

alphabet = 'NT'
figure()
t = ['Frequent mutations, N>', num2str(thres)];
plotFrequentMutations(frequentMutations, frequentMutationsCount, t, alphabet)

%% Analyze individual mutations
% mutations were previously found at
% nt1397 (nsp2), nt2891 (nsp3), nt3036 (=3037, nsp3), nt8782 (nsp4), nt11083 (nsp6),
% nt14408 (RdRp), nt17746 (=nt17747, helicase), nt17857 (=nt17858, helicase), 
% nt18060 (ExoN), nt23403 (S), nt26143 (ORF3a), nt28144 (ORF8a), nt28881 (N).
nt = 23403;                                                                 % Enter the nucleotide position of the mutation of interest
nt_protein = findProtein(nt);

j=0;
l=0;
mutated = [];
notMutated = []
for i=1:size(mutationArray,1)                                               % Loop through viruses
    mutationVector = mutationArray(i,:);                                    % Vector with location of mutations for virus i 
    if ismember(nt, mutationVector)==true                                   % Check if the position nt was mutated in virus i
        j = j+1;
        mutated(j) = length(mutationVector(~isnan(mutationVector)));        % If yes, store the number of mutations in 'mutated'
    else
        l = l+1;
        notMutated(l) = length(mutationVector(~isnan(mutationVector)));     % If not, store the number of mutations in 'notMutated'
    end
end

% Use ANOVA to check if there is a significant difference between the
% number of mutations in the groups mutated and notMutated:
maxMut = max(length(mutated), length(notMutated));                          % Maximal length of the two groups
notMutated = [notMutated,NaN(1,maxMut-length(notMutated))];                 % Make vectors the same size by adding NaNs to the shortest vector.
mutated = [mutated,NaN(1,maxMut-length(mutated))];

p = anova1([notMutated;mutated]',{'not mutated','mutated'},'off');          % ANOVA test
plotMutationAnalysis(nt, nt_protein, mutated, notMutated, p)

%% Analyze alignments per continent
% Define continents:
Europe = {'Europe','Czech Republic','Finland','France','Germany','Greece','Italy',...
    'Netherlands','Spain','Sweden','Turkey'};
Asia = {'Asia','China','Hong Kong','India','Iran','Israel','Japan','Kazakhstan',...
    'Nepal','Pakistan','South Korea','Sri Lanka','Taiwan'};
SouthAmerica = {'SouthAmerica','Brazil','Colombia','Peru','Puerto Rico'};
NorthAmerica = {'NorthAmerica','USA'};
Oceania = {'Oceania','Australia','Malaysia','Thailand','Viet Nam'};
Africa = {'Africa','South Africa'};

continents = {Africa,Asia,Europe,SouthAmerica,NorthAmerica,Oceania};

%% Store the mutation data per continent 

for i=1:length(continents)
    currCont = continents{i}{1};
    Nseqs.(currCont) = 0;                                                   % keep track of the number of sequences in that continent.
    mutsPerContinent.(currCont).Locations = NaN(length(seqs), maxNm);       % Store mutation location in an array in the structure 'mutsPerContinent'.
    mutsPerContinent.(currCont).Nm = [];                                    % Store the number of mutations in a vector.
end

for i=1:length(seqs)                                                        % Loop through sequences
    currCont = seqs(i).Continent;                                           % Get continent of current sequence
    if ~isempty(currCont)                                                   % If there is information about the continent ...
        Nseqs.(currCont) = Nseqs.(currCont)+1;                              % keep track of the number of sequences in that continent.
        currMutations = seqs(i).Mutations;                                  % get the mutation information of that sequence,
        mutsPerContinent.(currCont).Nm = [mutsPerContinent.(currCont).Nm, ...
            currMutations.Nm];                                              % append the number of mutations to the Nm vector in mutsPerContinent,
        mutationVector = currMutations.Locations;                           % get mutation locations,
        mutationVector = [mutationVector, ...                             	% make the vector the right length,
            NaN(1,maxNm-length(mutationVector))];   

        mutsPerContinent.(currCont).Locations(i,:) = mutationVector;     	% put vector in the mutation array of that continent.
    end
end
%% plot hist per continent
Nbins = [1,22,11,8,24,30];                                                  % Number of bins for each continent
% order: Africa, Asia, Europe, SouthAmerica, NorthAmerica, Oceania. 
figure()
for i=1:length(continents)
    currCont = continents{i}{1};                                            % Get current content (first entry in continents{i})
    subplot(3,2,i)
    hist(mutsPerContinent.(currCont).Nm, Nbins(i));                         % Histogram, number of bins = Nbins(i).
    title([currCont,', N=',num2str(Nseqs.(currCont))])
    xlim([0 60])
    xlabel('Number of mutations')
    ylabel('Count')
end
suptitle('SARS-CoV2 mutations across the world')

%% Plot the frequent mutations per continent
thresholds = [0,10,5,5,100,5];                                              % Frequent mutation threshold for each continent.
% order: Africa, Asia, Europe, NorthAmerica, Oceania. 
figure()
for i=1:length(continents)
    currCont = continents{i}{1};
    thres = thresholds(i);
    
    [freqMuts, freqMutsCount] = ...                                         % Find the locations of frequent mutations and their counts.
        evaluateFrequentMutations(mutsPerContinent.(currCont).Locations, thres);
    subplot(3,2,i)
    t = [currCont,' (Nseqs=',num2str(Nseqs.(currCont)), ') , N>',num2str(thres)];
    plotFrequentMutations(freqMuts, freqMutsCount, t, alphabet)             % Plot them in a barplot
end
suptitle('Frequent mutations of SARS-CoV2 in different locations')

%% Time analysis
periods = {'FEB','MARCH','APR1','APR2','MAY'};                              % Define time periods

for i=1:length(periods)
    currPeriod = periods{i};
    Nseqs.(currPeriod) = 0;                                                 % keep track of the number of sequences in that period.
    mutsPerPeriod.(currPeriod).Locations = NaN(length(seqs), maxNm);        % Store mutation location in an array in the structure 'mutsPerPeriod'.
    mutsPerPeriod.(currPeriod).Nm = [];                                     % Store the number of mutations in a vector.
end

for i=1:length(seqs)                                                        % Loop through sequences
    currPeriod = seqs(i).Period;                                            % Get period of current sequence
    if ~isempty(currPeriod)                                                 % If there is information about the continent ...
        Nseqs.(currPeriod) = Nseqs.(currPeriod)+1;                          % keep track of the number of sequences in that period.
        currMutations = seqs(i).Mutations;                                  % get the mutation information of that sequence,
        mutsPerPeriod.(currPeriod).Nm = [mutsPerPeriod.(currPeriod).Nm, ...
            currMutations.Nm];                                              % append the number of mutations to the Nm vector in mutsPerPeriod,
        mutationVector = currMutations.Locations;                           % get mutation locations,
        mutationVector = [mutationVector, ...                             	% make the vector the right length,
            NaN(1,maxNm-length(mutationVector))];   

        mutsPerPeriod.(currPeriod).Locations(i,:) = mutationVector;     	% put vector in the mutation array of that period.
    end
end

%% Plot hist per period
Nbins = [18,10,10,50,10];                                                    % Number of bins for each period.
% order: 'FEB','MARCH','APR1','APR2','MAY'. 
figure()
for i=1:length(periods)
    currPeriod = periods{i};
    subplot(3,2,i)
    hist(mutsPerPeriod.(currPeriod).Nm, Nbins(i));                          % plot histogram
    title([currPeriod,', N=',num2str(Nseqs.(currPeriod))])
    xlim([0 60])
    xlabel('Number of mutations')
    ylabel('Count')
end
suptitle('SARS-CoV2 mutations per time period')

%% Plot the frequent mutations per time period
thresholds = [1,1,50,50,20];                                                % Frequent mutation threshold for each time period.
% order: 'FEB','MARCH','APR1','APR2','MAY'. 
figure()
for i=1:length(periods)
    currPeriod = periods{i};
    thres = thresholds(i);
    
    [freqMuts, freqMutsCount] = ...                                         % Find the locations of frequent mutations and their counts.
        evaluateFrequentMutations(mutsPerPeriod.(currPeriod).Locations, thres);
    subplot(3,2,i)
    t = [currPeriod,' (Nseqs=',num2str(Nseqs.(currPeriod)), ') , N>',num2str(thres)];
    plotFrequentMutations(freqMuts, freqMutsCount, t, alphabet)             % Plot them in a barplot
end
suptitle('Frequent mutations of SARS-CoV2 in different time periods')

%% Protein sequence analysis
% To extract protein sequences and align them to the reference, see
% 'proteinAlignment.m'.
for s = 2:length(seqs)
    mismatches = char();                                                    % Keep track of the amino acid mismatches.
    proteinNames = seqs(s).Translation.Names;                               % Cell array with the names of the proteins of which a sequence is available for this virus.
    for p = 1:length(proteinNames)                                          % Loop through the available proteins
        alignment = seqs(s).Translation.Alignment{p};                       % Extract the alignment with the reference
        AA_Mismatches = char();                                             % Initialize mismatch array
        locations = [];                                                     % Ininialize a vector with locations of mismatches
        N_m = 0;                                                            % Keep track of the number of mutations
        N_mSimilar=0;                                                       % Keep track of the number of similar mutations (AAs in same group)
        N_mDifferent=0;                                                     % Keep track of the number of different mutations (AAs in different group)
        for i=1:size(alignment,2)                                           % Loop through alignment
            if alignment(1,i) ~= alignment(3,i) & alignment(1,i) ~= '-' ... % check if there is an AA difference 
                        & alignment(3,i) ~= 'X' & alignment(3,i) ~= '-'
                    if alignment(2,i) == ':'                                % Check if mismatch is similar (':' in alignment)
                        similarity = ' (similar)';
                        N_mSimilar = N_mSimilar+1;
                    else
                        similarity = '';
                        N_mDifferent = N_mDifferent+1;
                    end
                    N_m = N_m+1;
                    AA_Mismatches{N_m,1} = [num2str(i),': ', ...            % Store mismatches in format 'location: A->Y'
                        alignment(1,i),'->',alignment(3,i),similarity];
                    locations = [locations, i];                             % Store mismatch locations 
            end
        end
        mlocations{p} = locations;
        mismatches{p} = AA_Mismatches;
        Nmismatches(p) = N_m;
        NsimilarMismatches(p) = N_mSimilar;
        NdifferentMismatches(p) = N_mDifferent;
    end
    seqs(s).Translation.MismatchLocations = mlocations;                     % Store locations of mismatches
    seqs(s).Translation.Mismatches = mismatches;                            % Store mismatches
    seqs(s).Translation.Nm = Nmismatches;                                   % Store number of mismatches
    seqs(s).Translation.NmSimilar = NsimilarMismatches;                     % Store number of similar mismatches     
    seqs(s).Translation.NmDifferent = NdifferentMismatches;                 % Store number of different mismatches     
end

%% Histogram of proteins
% Define protein names of the reference.
proteinNamesRef = seqs(1).Translation.Names;
for i=1:length(proteinNamesRef)
    proteinNamesRefLower{i} = lower(proteinNamesRef{i});                    % List with lower case names
end
proteinAlternativeNamesRef = {'na','na','"spike glycoprotein"','na','na', ...	% Some proteins are named differently
    '"membrane protein"','na','na','na','na','"nucleocapsid protein"','na'};

% Initialize structure:
proteinAbbrevs = {'ORF1ab','ORF1a','Spike','ORF3a','Env','Membr','ORF6',...
    'ORF7a','ORF7b','ORF8','Nucl','ORF10'};
for p=1:length(proteinNamesRef)
    proteinAbbrev = proteinAbbrevs{p};
    misMatches.(proteinAbbrev) = [];                                        % Keep track of the number of mismatches in each protein
    sMisMatches.(proteinAbbrev) = [];                                       % Keep track of the number of similar mismatches in each protein
    dMisMatches.(proteinAbbrev) = [];                                       % Keep track of the number of different mismatches in each protein
end
% Find the number of mismatches
for s=1:length(seqs)
    proteinNames = seqs(s).Translation.Names;
    for p=1:length(proteinNames)
        protein = proteinNames{p};
        if (contains(protein,proteinNamesRef))==true | ...               	% If there is a corresponding protein in the reference           
                (contains(protein,proteinNamesRefLower))==true | ...        % (can also be in lower case letters), ...
                (contains(protein,proteinAlternativeNamesRef))==true        % (can also be one of the alternative names), ...
            m=1;
            while strcmp(protein, proteinNamesRef{m})==false & ...          % Find the index 'm' of the corresponding protein in the reference.
                    strcmp(protein, proteinNamesRefLower{m})==false & ...
                    strcmp(protein, proteinAlternativeNamesRef{m})==false
                m = m+1;
            end
            proteinAbbrev = proteinAbbrevs{m};
            misMatches.(proteinAbbrev) = ...                                % Store similar and non-similar mismatches
                [misMatches.(proteinAbbrev), seqs(s).Translation.Nm(p)];
            sMisMatches.(proteinAbbrev) = ...
                [sMisMatches.(proteinAbbrev), seqs(s).Translation.NmSimilar(p)];
            dMisMatches.(proteinAbbrev) = ...
                [dMisMatches.(proteinAbbrev), seqs(s).Translation.NmDifferent(p)];
        end
    end
end
%%
figure()
for p=1:length(proteinAbbrevs)
    proteinAbbrev = proteinAbbrevs{p};
    binrng = linspace(0,10,11);
    Nm = misMatches.(proteinAbbrev);
    
    subplot(3,4,p)
    histogram(Nm,binrng)
    xlim([0,10])
    xlabel('Number of nonsynonymous mutations')
    ylabel('Count')
    title([proteinAbbrev,', N=',num2str(length(Nm))]);
    hold off
end
suptitle('Nonsysnonymous mutations')

%% Evaluate frequent mismatches
% Find maximal number of mismatches
maxMisMatches = 0;
for p=1:length(proteinNamesRef)
    proteinAbbrev = proteinAbbrevs{p};
    if max(misMatches.(proteinAbbrev)) > maxMisMatches
        maxMisMatches = max(misMatches.(proteinAbbrev));
    end
end
% Initialize data structure mismatchesPerProtein with NaNs:
for p=1:length(proteinNamesRef)
    proteinAbbrev = proteinAbbrevs{p};
    mismatchesPerProtein.(proteinAbbrev) = NaN(length(seqs), maxMisMatches);
end
% Store the locations of mutations per protein: 
for s=1:length(seqs)
    proteinNames = seqs(s).Translation.Names;
    for p=1:length(proteinNames)
        protein = proteinNames{p};
        if (contains(protein,proteinNamesRef))==true | ...               	% If there is a corresponding protein in the reference           
                (contains(protein,proteinNamesRefLower))==true | ...        % (can also be in lower case letters), ...
                (contains(protein,proteinAlternativeNamesRef))==true        % (can also be one of the alternative names), ...
            m=1;
            while strcmp(protein, proteinNamesRef{m})==false & ...          % Find the index 'm' of the corresponding protein in the reference.
                    strcmp(protein, proteinNamesRefLower{m})==false & ...
                    strcmp(protein, proteinAlternativeNamesRef{m})==false
                m = m+1;
            end
        	proteinAbbrev = proteinAbbrevs{m};
            mismatchVector = seqs(s).Translation.MismatchLocations{p};    	% mismatchVector contains locations of mismatches
            mismatchVector = [mismatchVector, ...                           % Make mismatchVector the right length
                NaN(1, maxMisMatches-length(mismatchVector))];
            mismatchesPerProtein.(proteinAbbrev)(s,:) = mismatchVector;     % Store mismatchVector in mismatchesPerProtein.
        end
    end
end

alphabet = 'AA';                                                            % Amino acid alignment
thresholds = [100,50,5,50,0,1,0,1,0,10,10,0];                               % Frequent mismatch thresholds for all proteins
% Order: 'ORF1ab','ORF1a','Spike','ORF3a','Env','Membr','ORF6','ORF7a','ORF7b','ORF8','Nucl','ORF10'
for p=1:length(proteinNamesRef)                                             
    proteinAbbrev = proteinAbbrevs{p};  
    [freqMisMatches, freqMisMatchesCount] = ...                             % Evaluate frequent mismatches
            evaluateFrequentMutations(mismatchesPerProtein.(proteinAbbrev), thresholds(p));
    t = [proteinAbbrev,', N>',num2str(thresholds(p))];                      % Title of plot
    subplot(3,4,p)
    plotFrequentMutations(freqMisMatches, freqMisMatchesCount, t, alphabet) % Plot frequent mutations and their counts
end
suptitle('Frequent Amino Acid mismatches')

%% Functions

function uniqueMutations = findUniqueMutations(mutation_array)
% This function takes as input an array with the locations of mutations in
% all viruses.
% It outputs an array with 2 columns: the first column contains all unique
% mutation locations, the second column the frequency with which that
% location is mutated.
    [~,~,ic] = unique(mutation_array(:));                                       % Unique Values By Row, Retaining Original Order
    h = accumarray(ic, 1);                                                      % Count Occurrences
    counts = h(ic);                                                             % Map Occurrences To ‘ic’ Values
    uniqueMutations = unique([mutation_array(:),counts],'rows');
end

function protein = findProtein(nt)
% This function is a map of the SARS-CoV2 genome.
% It takes as input a nucleotide (integer between 1 and 29903).
% It outputs in which open reading frame the nucleotide is present.
    if nt>=1 && nt<=265
        protein = '5-UTR';
    elseif nt>=266 && nt<=805
        protein = 'nsp1';
    elseif nt>=806 && nt<=2719 
        protein = 'nsp2';
    elseif nt>=2720 && nt<=8554
        protein = 'nsp3';
    elseif nt>=8555 && nt<=10054
        protein = 'nsp4';
    elseif nt>=10055 && nt<=10972
        protein = 'nsp5';
    elseif nt>=10973 && nt<=11842
        protein = 'nsp6';
    elseif nt>=11843 && nt<=12091
        protein = 'nsp7';
    elseif nt>=12092 && nt<=12685
        protein = 'nsp8';
    elseif nt>=12686 && nt<=13024
        protein = 'nsp9';
    elseif nt>=13025 && nt<=13441
        protein = 'nsp10';
    elseif nt>=13442 && nt<=16236
        protein = 'RdRp';
    elseif nt>=16237 && nt<=18039
        protein = 'helicase';
    elseif nt>=18040 && nt<=19620
        protein = 'ExoN';
    elseif nt>=19621 && nt<=20658
        protein = 'endoRNAse';
    elseif nt>=20659 && nt<=21552
        protein = 'methyltransferase';
    elseif nt>=21563 && nt<=25384
        protein = 'S';
    elseif nt>=25393 && nt<=26220
        protein = 'ORF3a';
    elseif nt>=26245 && nt<=26472
        protein = 'E';
    elseif nt>=26523 && nt<=27191
        protein = 'M';
    elseif nt>=27202 && nt<=27387
        protein = 'ORF6';
    elseif nt>=27394 && nt<=27759
        protein = 'ORF7a';
    elseif nt>=27756 && nt<=27887
        protein = 'ORF7b';
    elseif nt>=27894 && nt<=28259
        protein = 'ORF8';
    elseif nt>=28274 && nt<=29533
        protein = 'N';
    elseif nt>=29558 && nt<=29674
        protein = 'ORF10';
    elseif nt>=29675 && nt<=29903
        protein = '3-UTR';
    else
        protein = 'None';
    end
end

function plotFrequentMutations(frequentMutations, countFrequentMutations, t, alphabet)
% This function plots a barplot of the most frequent mutations. t=title
    for i=1:length(frequentMutations)                                           % Loop through frequend mutations
        mutated_nt = frequentMutations(i);                                      
        protein = findProtein(mutated_nt);                                      % Find in which ORF the mutation is present
        if strcmp(alphabet, 'NT')
            x_label{i} = join([num2str(mutated_nt),' (',protein,')']);       	% Store it in a cell array
            xAxis = 'Location in the genome';
        else
            x_label{i} = num2str(mutated_nt);
            xAxis = 'Location in protein';
        end
    end
    N_frequentMutations = length(countFrequentMutations);
    bar(1:N_frequentMutations, countFrequentMutations)
    set(gca,'XTickLabel',x_label)
    set(gca,'XTick',1:N_frequentMutations)
    set(gca,'Fontsize',12)
    set(gca,'TickLabelInterpreter', 'tex');
    xtickangle(45)
    title(t)
    xlabel(xAxis)
    ylabel('Count')
    hold off
end

function [frequentMutations, frequentMutationsCount] = evaluateFrequentMutations(mutationArray, thres)
    uniqueMutations = findUniqueMutations(mutationArray);                       % Find which unique mutations are present and with what frequencies
    allMutations = uniqueMutations(:,1);                                        % Unique mutations
    allMutationsCount = uniqueMutations(:,2);                                   % Their frequencies
    allMutations = allMutations(~isnan(allMutations));
    allMutationsCount = allMutationsCount(~isnan(allMutations));

    frequentMutations = allMutations(allMutationsCount>thres);                  % Filter mutations which are less frequent than thres
    frequentMutationsCount = allMutationsCount(allMutationsCount>thres);        % Find their frequencies
end

function plotMutationAnalysis(nt, protein, mutated, notMutated, p)
% This function plots the result of the analysis of individual mutations.
    figure()
    suptitle(['Analysis of mutation nt',num2str(nt),' (',protein,')'])
    N_mut = length(mutated(~isnan(mutated)));
    N_notMut = length(notMutated(~isnan(notMutated)));
    subplot(1,2,1)
    histogram(notMutated)
    hold on
    histogram(mutated)
    xlabel('Number of mutations')
    ylabel('Count')
    
    l1 = {['nt',num2str(nt),' not mutated, N=',num2str(N_notMut)],...
        ['nt',num2str(nt),' mutated, N=',num2str(N_mut)]};
    legend(l1)
    title('Histogram')
    hold off

    subplot(1,2,2)
    l2 = {['nt',num2str(nt),' not mutated'],['nt',num2str(nt),' mutated']};
    boxplot([notMutated;mutated]','Labels',l2)
    xticks = [1,2];
    ylabel('Number of mutations')
    title(['Boxplot, p=',num2str(p)])
    hold off
end
