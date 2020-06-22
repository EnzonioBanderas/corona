function [seq_loc, seq_mut, seq_nMut, ...
    aseq_loc, aseq_mut, aseq_nMut, ...
    aseq_alltime_unique_n, ...
    V, I, r, d, dtot, ntot, nAA] = ...
remove(i0, ...
    seq_loc, seq_mut, seq_nMut, ...
    aseq_loc, aseq_mut, aseq_nMut, ...
    aseq_alltime_unique_loc, aseq_alltime_unique_mut, aseq_alltime_unique_n, ...
    V, I, r, d, dtot, ntot, nAA)

    oldAAsequence_loc = aseq_loc{i0};
    oldAAsequence_mut = aseq_mut{i0};
    oldAAsequence_nMut = aseq_nMut(i0);
    seq_loc{i0} = []; seq_mut{i0} = [];
    aseq_loc{i0} = []; aseq_mut{i0} = []; 
    seq_nMut(i0) = 0; aseq_nMut(i0) = 0;
    for iAA = find(aseq_alltime_unique_n~=0)
        current_loc = aseq_alltime_unique_loc{iAA};
        current_mut = aseq_alltime_unique_mut{iAA};
        if oldAAsequence_nMut==length(current_loc)
            if all(oldAAsequence_loc==current_loc)
                if all(oldAAsequence_mut==current_mut)
                    aseq_alltime_unique_n(iAA) = aseq_alltime_unique_n(iAA) - 1;
                    if aseq_alltime_unique_n(iAA)==0
                        nAA = nAA - 1;
                    end
                    break
                end   
            end
        end
    end
    r(i0) = 0;
    d(:,i0) = zeros(size(d,1), 1);
    dtot(i0) = 0;
    ntot = ntot - 1;
end              