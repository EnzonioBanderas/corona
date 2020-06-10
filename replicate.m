function [seq, aseq, ntot, V, r, I, d, dtot] = ...
    replicate(i0, seq, aseq, mu, ntot, V, r, I, d, dtot, sigma, r0, refSeq, pInfo)
% This function lets a virus of strain i0 replicate. During the
% replication, each nucleotide has a probability mu to mutate.
% If a new strain is formed by mutation, the fitness of this new strain is
% calculated with the fitness function.
    oldsequence = seq{i0};
    newsequence = oldsequence;
    L = length(oldsequence);
    nt = ['a','t','g','c'];
    % loop through nucleotides in sequence
    mutations = rand(1,L);
    mutationBoolean = (mutations <= mu);
    for n = find(mutationBoolean)
        % determine if a mutation occurs at this location:
      	newnt = nt(nt ~= newsequence(n));
       	newsequence(n) = newnt(randi(3));
    end
    if strcmp(oldsequence, newsequence)
        V(i0) = V(i0) + 1;
        return
    else
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
            aseq{newLoc} = nt2aa(newsequence);
            V(newLoc) = 1;
            di = distance(refSeq, aseq{newLoc}, pInfo);
            d(:,newLoc) = di;
            dtot(newLoc) = sum(di);
            r(newLoc) = fitness(di, r0, sigma, pInfo);
            I(newLoc) = 0; 
        end
    end
end