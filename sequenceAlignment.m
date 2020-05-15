%{
12-04-2020

This matlab file:
- finds SARS-CoV2 GenBank sequences which have between 28000 and 30000 resolved nucleotides;
- loads these (2277) sequences from GenBank, takes ~1.5 hrs;
- aligns all these sequences to the reference sequence (NC_045512 in genbank), takes ~5 hrs;
- saves the sequence ID's (names) and data (seqs) to Sequence_alignment.mat.

TOTAL RUNTIME: ~6.5 hrs
%}

%% Load data
clear all
close all

% file = fasta file with all 2434 full SARS-CoV2 sequences. 
% Can be downloaded from 'https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049'
file = 'sequences.fasta.txt';
data = fastaread(file);                                    

%% Find sequences that have between 28000 and 30000 resolved nucleotides:
j=0;
for i = 1:length(data)
    L = length(data(i).Sequence);                                           % length of sequence
    N = count(data(i).Sequence, 'N');                                       % number of unresolved nucleotides ('N' in fasta file)
    if L-N > 28000 & L-N < 30000                                            % check if the number of unresolved nucleotides (L-N) is between 28000 and 30000.
        j = j+1;
        header = split(data(i).Header, ' |');
        names(j) = header(1);                                               % store the name / id of the genbank entry
    end
end

%% Load sequences in 'names' from GenBank
for i = 1:length(names)
    seqs(i) = getgenbank(names{i});
    if mod(i,100) == 0
        disp(i)
    end
end

%% Needleman-Wunsch Glocal alignment
% In a "glocal" alignment, gap penalties at the end of the sequences are null.
% Aligning 2 sequences takes ~6 seconds.
Seq1 = seqs(1).Sequence;                                                    % reference sequence (NC_045512 in genbank)
for i=1:length(seqs)
    Seq2 = seqs(i).Sequence;                                                % other sequence 
    [score, alignment, start] = ...                                         % needleman-Wunsch alignment
        nwalign(Seq1,Seq2, 'Alphabet', 'NT','Glocal',true);
    seqs(i).('Alignment') = alignment;                                      % Store sequence alignment in structure 'seqs'
    seqs(i).('Score') = score;                                              % Needleman-Wunsch alignment score.
    if mod(i,100) == 0
        disp(i)
    end
end

%% Save results
save('Sequence_alignment.mat','names','seqs')
