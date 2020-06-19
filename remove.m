function [seq_loc, seq_mut, aseq_loc, aseq_mut, V, I, r, d, dtot, ntot, nAA, seq_nMut] = ...
             remove(i0, seq_loc, seq_mut, aseq_loc, aseq_mut, V, I, r, d, dtot, ntot, nAA, seq_nMut)

%     oldAAsequence = aseq{i0};
    seq_loc{i0} = []; seq_mut{i0} = [];
    aseq_loc{i0} = []; aseq_mut{i0} = []; 
    seq_nMut(i0) = 0;
%     sameSeq = find(strcmp(aseq, oldAAsequence),1);
%     if isempty(sameSeq)     % the AA sequence is no longer present
%         nAA = nAA - 1;
%     end 
    r(i0) = 0;
    d(:,i0) = zeros(size(d,1), 1);
    dtot(i0) = 0;
    ntot = ntot - 1;
end              