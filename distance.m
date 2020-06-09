function dist = distance(refseq, seq, pInfo)
% This function calculates the Hamming distance between the reference
% sequence and another sequence, for all proteins.
% dist = a column vector with distances of all proteins.
    
    pNames = fields(pInfo);
    dist = zeros(length(pNames),1);
    
    for i = 1:length(pNames)
        protein = pNames{i};
        startP = pInfo.(protein).proteinLocation(1);
        endP = pInfo.(protein).proteinLocation(2);
                
        pRefSeq = refseq(startP:endP);
        pSeq = seq(startP:endP);
        
        dist(i) = sum(pRefSeq ~= pSeq);
    end
end