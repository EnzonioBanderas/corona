function [seq, aseq, ntot, V, r, I, d] = ...
    replicate(i0, seq, aseq, mu, ntot, V, r, I, d, sigma, r0, refSeq, pInfo)
% This function lets a virus of strain i0 replicate. During the
% replication, each nucleotide has a probability mu to mutate.
% If a new strain is formed by mutation, the fitness of this new strain is
% calculated with the fitness function.
    oldsequence = seq{i0};
    newsequence = oldsequence;
    L = length(oldsequence);
    nt = ['a','t','g','c'];
    % loop through nucleotides in sequence
    for n = 1:L
        % determine if a mutation occurs at this location:
        if rand <= mu
            newnt = nt(nt ~= newsequence(n));
            newsequence(n) = newnt(randi(3));
        end
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
            seq{ntot} = newsequence;
            aseq{ntot} = nt2aa(newsequence);
            V(ntot) = 1;
            di = distance(refSeq, aseq{ntot}, pInfo);
            d(:,ntot) = di;
            r(ntot) = fitness(di, r0, sigma, pInfo);
            I(ntot) = 0;
        end
    end
end