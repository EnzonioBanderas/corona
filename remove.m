function [seq, aseq, V, I, r, d, dtot, ntot] = ...
             remove(i0, seq, aseq, V, I, r, d, dtot, ntot)
    %extinct = (V+I == 0);  
    %if any(extinct)
%     newseq = cell(size(seq));
%     newseq(~extinct) = seq(~extinct);
% 
%     newaseq = cell(size(aseq));
%     newaseq(~extinct) = aseq(~extinct);
    seq{i0} = {};
    aseq{i0} = {}; 
    r(i0) = 0;
    d(:,i0) = zeros(size(d,1), 1);
    dtot(i0) = 0;
    ntot = ntot - 1;
    %end
end
