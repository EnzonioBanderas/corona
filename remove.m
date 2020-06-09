function [seq, aseq, V, I, r, d, ntot] = ...
             remove(seq, aseq, V, I, r, d, ntot)
    extinct = (V+I == 0);  
    if any(extinct)
        V(extinct) = [];
        I(extinct) = [];
        seq(extinct) = [];
        aseq(extinct) =[];
        r(extinct) = [];
        d(:,extinct) = [];
        ntot = ntot - sum(extinct);
    end
end
