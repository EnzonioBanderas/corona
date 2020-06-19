function [seq, aseq, ntot, nAA, V, r, I, d, dtot] = ...
    replicate(i0, seq, aseq, mu, ntot, nAA, V, r, I, d, dtot, sigma, r0, refSeq, pInfo)
% This function lets a virus of strain i0 replicate. During the
% replication, each nucleotide has a probability mu to mutate.
% If a new strain is formed by mutation, the fitness of this new strain is
% calculated with the fitness function.
    oldsequence = seq{i0};
    newsequence = oldsequence;
    L = length(oldsequence);
    nt = ['A','T','G','C'];
    
    mutations = rand(1,L);
    mutationBoolean = (mutations <= mu);
    % Loop through mutation positions:
    for n = find(mutationBoolean)
      	newnt = nt(nt ~= newsequence(n));
       	newsequence(n) = newnt(randi(3));
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
        newLoc = find(V + I == 0, 1);                                   % location of new strain
        seq{newLoc} = newsequence;

        oldAAsequence = aseq{i0};
        newAAsequence = oldAAsequence;
        % translate the protein
        for n = find(mutationBoolean)
            if mod(n, 3) == 0       % we are at the end of a codon
                aLoc = n / 3;       % protein position 
                codon = newsequence(n-2:n);
            elseif mod(n, 3) == 1   % we are at the beginning of a codon
                aLoc = ceil(n / 3); % protein position 
                codon = newsequence(n:n+2);
            elseif mod(n, 3) == 2   % we are in the middle of a codon
                aLoc = ceil(n / 3); % protein position
                codon = newsequence(n-1:n+1);
            end
            newAA = nt2aa(codon);
            newAAsequence(aLoc) = newAA;
        end
%         sameSeq = find(strcmp(aseq, newAAsequence),1);
%         if isempty(sameSeq)
%             nAA = nAA + 1;
%         end
        aseq{newLoc} = newAAsequence;
        V(newLoc) = 1;
        I(newLoc) = 0; 
        di = distance(refSeq, aseq{newLoc}, pInfo);
        d(:,newLoc) = di;
        dtot(newLoc) = sum(di);
        r_i0 = r(i0);
        if strcmp(newAAsequence, oldAAsequence)
            r(newLoc) = r_i0;
        else
            r(newLoc) = replicationRate(di, r_i0, sigma, pInfo);
        end   
    end
    
end