function [seq, aseq, V, I, r, d, dtot, ntot, nAA] = ...
             remove(i0, seq, aseq, V, I, r, d, dtot, ntot, nAA)

    oldAAsequence = aseq{i0};
    seq{i0} = {};
    aseq{i0} = {}; 
%     sameSeq = find(strcmp(aseq, oldAAsequence),1);
%     if isempty(sameSeq)     % the AA sequence is no longer present
%         nAA = nAA - 1;
%     end 
    r(i0) = 0;
    d(:,i0) = zeros(size(d,1), 1);
    dtot(i0) = 0;
    ntot = ntot - 1;
end
