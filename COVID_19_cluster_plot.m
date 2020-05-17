clearvars

load('Sequence_alignment.mat')
seqs=seqs(1:200);
nSeqs = length(seqs);

mutation_table_cell = cell(1, nSeqs);
for iSeqs=1:nSeqs
    
%     iSeqs
%     nAlign = size(seqs(iSeqs).Alignment,2);
    iAlign_index = find((seqs(iSeqs).Alignment(1,:)~=seqs(iSeqs).Alignment(3,:)));
%     &...
%         ((seqs(iSeqs).Alignment(1,:)~='-')&(seqs(iSeqs).Alignment(3,:)~='-'&seqs(iSeqs).Alignment(3,:)~='N')));
    nAlign_index = length(iAlign_index);
    mutation_table_ref = cell(nAlign_index, 1);
    mutation_table_ind = zeros(nAlign_index, 1);
    mutation_table_seq = mutation_table_ref;
    for iAlign=1:nAlign_index
        mutation_table_ref{iAlign} = seqs(iSeqs).Alignment(1, iAlign_index(iAlign));
        mutation_table_ind(iAlign) = iAlign_index(iAlign);
        mutation_table_seq{iAlign} = seqs(iSeqs).Alignment(3, iAlign_index(iAlign));
    end
    mutation_table_cell{iSeqs} = table(mutation_table_ref,...
        mutation_table_ind, mutation_table_seq,...
        'VariableNames', {'ref', 'ind', 'seq'});
    
end

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

tic
nRow = length(mutation_table_cell);
seqs_distances = zeros(nRow);
for i=1:nRow
    for j=1:i
%         i
%         j
        [Mu,ia,ic] = unique([mutation_table_cell{i}; mutation_table_cell{j}], 'rows', 'stable'); % Unique Values By Row, Retaining Original Order
        h = accumarray(ic, 1); % Count Occurrences
        maph = h(ic); % Map Occurrences To ‘ic’ Values
        seqs_distances(i, j) = sum(maph==1);
        seqs_distances(j, i) = seqs_distances(i, j);
    end
end
toc

tic
coord_MDS = cmdscale(seqs_distances, 2);
toc

scatter(coord_MDS(:, 1), coord_MDS(:, 2), 10, categorical({seqs.Country}), 'filled')
title('Sequence scatter plot, dimensional reduction with MDS on Manhattan distance matrix')
xlabel('MDS 1')
ylabel('MDS 2')