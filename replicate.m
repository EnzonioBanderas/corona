function [seq_loc, seq_mut, aseq_loc, aseq_mut, ntot, nAA, V, r, I, d, dtot, seq_nMut] = ...
    replicate(myStream, i0, seq_loc, seq_mut, aseq_loc, aseq_mut, mu, ntot, nAA, V, r, I, d, dtot, sigma, r0, gRefSeq, L, pRefSeq, pInfo, seq_nMut, proteinLocation, translateCodon)
% This function lets a virus of strain i0 replicate. During the
% replication, each nucleotide has a probability mu to mutate.
% If a new strain is formed by mutation, the fitness of this new strain is
% calculated with the fitness function.
    oldsequence_loc = seq_loc{i0};
    oldsequence_mut = seq_mut{i0};
    nt = ['A','T','G','C'];
    
    % Get mutation locations and number of mutations
    mutations = rand(myStream,1,L);
    mutationBoolean = (mutations <= mu);
    loc = find(mutationBoolean); % find output is always sorted!
    nMut = length(loc);
    
    % If new sequence is old sequence, add viral particle and exit
    if nMut==0 % newsequence is oldsequence
        V(i0) = V(i0) + 1;
        return
    end
    
    % Mutate new sequence by combining old sequence mutations with new
    % sequence mutations.
    iMutRandi = randi(myStream,3, [1, nMut]);
    mut = blanks(nMut);
    % Go through mutation locations, if it was mutated in old mutations
    % then use the old mutation to get the 3 possible nucleotides to which
    % the position can mutate, else if it is a previously unmutated
    % position get the nucleotide of the reference to determine the 3
    % possible nucleotides.
    newpos_in_oldpos = ismember(loc, oldsequence_loc);
    for iMut = 1:nMut
        if newpos_in_oldpos(iMut) % position was mutated in old sequence, use old mutations
            newnt = nt(nt ~= oldsequence_mut(oldsequence_loc == loc(iMut)));
        else % position was not mutated in old sequence, use reference
            newnt = nt(nt ~= gRefSeq(loc(iMut)));
        end
       	mut(iMut) = newnt(iMutRandi(iMut));
    end
    % Merge old mutations with new mutations
    oldpos_not_in_newpos = ~ismember(oldsequence_loc, loc);
    newsequence_loc = [oldsequence_loc(oldpos_not_in_newpos), loc]; % remove oldsequence locs which are in loc
    newsequence_mut = [oldsequence_mut(oldpos_not_in_newpos), mut]; % remove oldsequence locs which are in loc

    % Sort after merging old and new mutations
    [newsequence_loc, sortmut] = sort(newsequence_loc); 
    newsequence_mut = newsequence_mut(sortmut);
    
    % Correct newsequence for backmutation to reference 
    backmutation_correction = gRefSeq(newsequence_loc)~=newsequence_mut;
    newsequence_loc = newsequence_loc(backmutation_correction);
    newsequence_mut = newsequence_mut(backmutation_correction);
    
    % Get length of newsequence, if 0 then add viral particle to reference
    nMut = length(newsequence_loc);
    seq_nMut0 = seq_nMut == 0; % empty sequences without mutations, first one is reference
    if nMut==0
        iRef = find(seq_nMut0, 1);
        V(iRef) = V(iRef) + 1;
        return
    end
    
    % Go through all sequences and check whether locations are shared
    % If locations are shared, check whether mutations are shared
    % If both of these arrays are exactly the same, add viral particle to
    % the current sequence in the loop and return
    seq_nMutNot0_index = find(~seq_nMut0);
    for iSeq = seq_nMutNot0_index
        % If number of mutations is the same
        if seq_nMut(iSeq)==nMut
            % if all positions are the same between old and new
            if all(seq_loc{iSeq}==newsequence_loc)
                % if also all mutations are the same between old and new, add viral
                % particle to current sequence in loop (sorting is important)
                if all(seq_mut{iSeq}==newsequence_mut)
                    V(iSeq) = V(iSeq) + 1;
                    return
                end
            end
        end
    end

    % newsequence is truly new, add it to sequence cell array
    ntot = ntot + 1;
    newLoc = find(V + I == 0, 1); % location of new strain
    seq_loc{newLoc} = newsequence_loc;
    seq_mut{newLoc} = newsequence_mut;
    seq_nMut(newLoc) = nMut;

    % Translate newsequence to AA sequence by looking at which codons are
    % affected in new mutations (info not in newsequence_loc but in loc).
    % Get unique array of codons which were newly affected and add old
    % mutations also affecting these newly affected codons, then loop over the unique codons
    % and take the reference codon for each newly affected codon and mutate
    % it according to the merged mutations of oldsequence and newsequence.
    % Then translate the mutated reference codon. Record the codon position and what amino
    % acid this position is mutated to.
    newly_affected_codon_unique = unique(ceil(loc/3)); % newly affected codon index (from oldsequence)
    codonIndex_seq = ceil(newsequence_loc/3); % mutation in newsequence affects this codon
    affected_codon_nt_index = ismember(codonIndex_seq, newly_affected_codon_unique);
    codonIndex_seq_affected = codonIndex_seq(affected_codon_nt_index);
    newsequence_loc_affected = newsequence_loc(affected_codon_nt_index);
    newsequence_mut_affected = newsequence_mut(affected_codon_nt_index);
    nCodon = length(newly_affected_codon_unique);
    mut_AA = blanks(nCodon);
    for iCodon = 1:nCodon % loop through newly affected codons
        % apply mutations on gRefSeq codons
        gRefSeqEndCodon = newly_affected_codon_unique(iCodon)*3;
        gRefSeqStartCodon = gRefSeqEndCodon - 2;
        gRefSeqCodon = gRefSeq(gRefSeqStartCodon:gRefSeqEndCodon);
        iCodon_affected_nt = codonIndex_seq_affected == newly_affected_codon_unique(iCodon); % select all mutations which affect current codon
        newsequence_loc_affected_iCodon = mod(newsequence_loc_affected(iCodon_affected_nt)-1, 3)+1;
        newsequence_mut_affected_iCodon = newsequence_mut_affected(iCodon_affected_nt);
        gRefSeqCodon(newsequence_loc_affected_iCodon) = newsequence_mut_affected_iCodon; % ref codon is mutated
        mut_AA(iCodon) = translateCodon.(gRefSeqCodon); % ref codon is translated (any benefit to vectorize this? gRefSeqCodon would become cell array)
    end
    oldAAsequence_loc = aseq_loc{i0};
    oldAAsequence_mut = aseq_mut{i0};
    newAAsequence_loc = [oldAAsequence_loc, newly_affected_codon_unique];
    newAAsequence_mut = [oldAAsequence_mut, mut_AA];
    
    % sort after merging old and new mutations (AA)
    [newAAsequence_loc, sortmut] = sort(newAAsequence_loc); 
    newAAsequence_mut = newAAsequence_mut(sortmut);
    
    % Correct newAAsequence for backmutation to reference 
    backmutation_correction_AA = pRefSeq(newAAsequence_loc)~=newAAsequence_mut;
    newAAsequence_loc = newAAsequence_loc(backmutation_correction_AA);
    newAAsequence_mut = newAAsequence_mut(backmutation_correction_AA);
    
    aseq_loc{newLoc} = newAAsequence_loc;
    aseq_mut{newLoc} = newAAsequence_mut;
    V(newLoc) = 1;
    I(newLoc) = 0;
    di = distance(newAAsequence_loc, proteinLocation);
    d(:,newLoc) = di;
    dtot(newLoc) = sum(di);
    r_i0 = r(i0);
    
    % Check whether amino acid sequence is the same as old
    if length(newAAsequence_loc)==length(oldAAsequence_loc)
        if all(newAAsequence_loc==oldAAsequence_loc)
            if all(newAAsequence_mut==oldAAsequence_mut)
                r(newLoc) = r_i0;
                return
            end   
        end
    end
    r(newLoc) = replicationRate(di, r_i0, sigma, pInfo);
    
%         sameSeq = find(strcmp(aseq, newAAsequence),1);
%         if isempty(sameSeq)
%             nAA = nAA + 1;
%         end
    
end