function [gRefSeq, pRefSeq, pNames, proteinLocation, genomeLocation] = getRefSeq()

    % Load data
    pNames = ...
        {'NSP1', 'NSP2', 'NSP3', 'NSP4', 'NSP5', 'NSP6', 'NSP7', 'NSP8', ...
        'NSP9', 'NSP10', 'NSP12', 'NSP13', 'NSP14', 'NSP15', 'NSP16', ...
        'Spike', 'NS3', 'E', 'M', 'NS6', 'NS7a', 'NS7b', 'NS8', 'N'};
    N = length(pNames);                                                     % number of proteins
    gData = load('Data/refSeqGenome.mat');
    gData = gData.refseq;
    
    gRefCell = genome2cell(gData);                                          % cell array with all genome sequences     
    pRefCell = cell(1, N);
    for p = 1:N
        pRefCell{p} = nt2aa(gRefCell{p}); 
    end
    % combine the cell arrays into two character arrays
    gRefSeq = [];                                                           % char array with all genome sequences 
    pRefSeq = [];                                                           % char array with all protein sequences
    for p = 1:length(gRefCell)
        pRefSeq = [pRefSeq, pRefCell{p}];
        gRefSeq = [gRefSeq, gRefCell{p}];
    end
    % Collect starting and ending locations of each protein
    proteinLocation = zeros(length(pNames), 2);
    genomeLocation = proteinLocation;
    startP = 1;                                                             % start of protein in protein sequence
    startG = 1;                                                             % start of protein in genome sequence
    for i = 1:length(pNames)
        L = length(pRefCell{i});                                            % length of protein
        proteinLocation(i, :) = [startP, min(startP + L -1, length(pRefSeq))];  % start and end of protein in amino acid sequence
        genomeLocation(i, :) = [startG, min(startG + 3*L, length(gRefSeq))];    % start and end of protein in nuclotide sequence
        startP = startP + L;
        startG = startG + 3*L;
    end
end