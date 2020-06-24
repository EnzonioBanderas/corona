function [seq_loc, seq_mut, seq_nMut, ...
    aseq_loc, aseq_mut, aseq_nMut, ...
    V, I, r, d, dtot, ntot, nAA] = ...
remove(i0, ...
    seq_loc, seq_mut, seq_nMut, ...
    aseq_loc, aseq_mut, aseq_nMut, ...
    V, I, r, d, dtot, ntot, nAA)

    oldAAsequence_loc = aseq_loc{i0};
    oldAAsequence_mut = aseq_mut{i0};
    oldAAsequence_nMut = aseq_nMut(i0);
    seq_loc{i0} = []; seq_mut{i0} = [];
    aseq_loc{i0} = []; aseq_mut{i0} = []; 
    seq_nMut(i0) = 0; aseq_nMut(i0) = 0;
    
    r(i0) = 0;
    d(:,i0) = zeros(size(d,1), 1);
    dtot(i0) = 0;
    ntot = ntot - 1;
    
    % Look through all remaining AAs, if the same AA is still present exit,
    % the function and do not decrease AA.
    remaining = find((V+I)~=0);
    for iAA = remaining 
        if oldAAsequence_nMut==aseq_nMut(iAA)
            if all(oldAAsequence_loc==aseq_loc{iAA})
                if all(oldAAsequence_mut==aseq_mut{iAA})
                    % The sequence was identical to another already
                    % existing sequence, so return as nAA remains the
                    % same.
                    return
                end   
            end
        end
    end
    nAA = nAA - 1; % AA not the same as any other AA, reduce nAA
    
end              