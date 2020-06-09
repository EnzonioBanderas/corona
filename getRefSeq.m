function [gRefSeq, pRefSeq, pNames, pInfo] = getRefSeq()

    % Load data
    pData = load('refSeqProtein.mat');
    gData = load('refSeqGenome.mat');
    gData = gData.refseq;
    pData = pData.proteinref.Translation;
    
    pNames = fields(pData);                                                 % array with names of proteins
    pNames(11) = [];                                                        % Remove NSP11 (is already in NSP12)
    pRefCell = cell(1, length(pNames));                                     % cell array with protein sequences
    % Fill in the reference protein sequence cell
    for i=1:length(pNames)
        protein = pNames{i};
        pRefCell{i} = pData.(protein).Sequence;
    end
    pRefCell{11}(9) = 'K';                                                  % Fix one substitution error                                                
    gRefCell = genome2cell(gData);                                          % cell array with all genome sequences                          
    % combine the cell arrays into two character arrays
    gRefSeq = [];                                                           % char array with all genome sequences 
    pRefSeq = [];                                                           % char array with all protein sequences
    for p = 1:length(gRefCell)
        pRefSeq = [pRefSeq, pRefCell{p}];
        gRefSeq = [gRefSeq, gRefCell{p}];
    end
    % Collect starting and ending locations of each protein
    pInfo = struct();                                                       % structure containing start and end of all proteins
    startP = 1;                                                             % start of protein in protein sequence
    startG = 1;                                                             % start of protein in genome sequence
    for i = 1:length(pNames)
        protein = pNames{i};
        L = length(pRefCell{i});
        pInfo.(protein).('proteinLocation') = [startP, min(startP + L -1, length(pRefSeq))];
        pInfo.(protein).('genomeLocation') = [startG, min(startG + 3*L, length(gRefSeq))];
        startP = startP + L;
        startG = startG + 3*L;
    end
end