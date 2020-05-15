%{
Extract protein sequences from 'seqs.Features' and perform a
needleman-Wunsch alignment with the reference sequence.
%}

%% Extract the protein sequences
for s = 1:length(seqs)
    info = seqs(s).Features;                                                % seqs(s).Features contains information about the proteins.
    j=0;
    k=0;
    proteinSequences = char();                                          	% Initialize array with protein sequences.
    proteinNames = char();                                                  % Initialize array with protein names.
    proteinSequence = false;                                                % Boolean to indicate whether a line in 'info' contains a protein sequence. Sequences usually span multiple lines.
    for i=1:size(info,1)                                                    % Loop through line in info.
        splitInfo = strsplit(info(i,:),'  ');                               % Split the line
        if proteinSequence == true                                          % If the line contains a protein sequence:
            splitString = split(splitInfo{2},'"');                          % Split the second part of the line ('"' always appears at the end of a sequence).
            proteinSequences{k} = [proteinSequences{k},splitString{1}];     % The first part of splitString contains the sequence. Append it to proteinSequences{k}.
            if length(splitString)>1                                        % If you are at the end of a sequence (then the line contained '"'), ...
                proteinSequence = false;                                    % the following line in info is no longer a protein sequence.
            end
        end
        if contains(splitInfo{2},'/translation=')                           % If the line in info starts with /translation=, we are at the beginning of a sequence.
            k = k+1;
            proteinSequence = true;                                         % From now on we are reading a protein sequence.
            splitString = split(splitInfo{2},'="');                         % Split the line on =". The second part contains the protein sequence.
            AAseq = split(splitString{2}, '"');                             % If the sequence is short and spans only 1 line, a " will appear at the end, ...
            if length(AAseq)>1                                              % and the next line will no longer be a sequence.
                proteinSequence = false;
            end
            proteinSequences{k} = AAseq{1};                                 % Start the kth sequence.
            % Extract protein name:
            c = 1;
            nameInfo = splitInfo;
            while contains(nameInfo{2},'product') == false                  % The protein name is a few lines above the translation line. It is indicated with the word 'product'.
                nameInfo = strsplit(info(i-c,:),'  ');
                c = c+1;
            end
            nameString = split(nameInfo(2),'=');
            proteinNames{k} = nameString{2};                                % Append the protein name to proteinNames.
        end
    end
    seqs(s).('Translation').('Sequences') = proteinSequences;               % Store protein names and sequences in seqs.Translation.
    seqs(s).('Translation').('Names') = proteinNames;
end

%% Show a histogram of how many proteins are present
Nproteins = [];
for i = 1:length(seqs)
    Nproteins = [Nproteins,length(seqs(i).Translation.Sequences)];
end
histogram(Nproteins)
%% Define protein names and sequences of the reference.
proteinSequencesRef = seqs(1).Translation.Sequences;
proteinNamesRef = seqs(1).Translation.Names;

for i=1:length(proteinNamesRef)
    proteinNamesRefLower{i} = lower(proteinNamesRef{i});                    % List with lower case names
end
proteinAlternativeNamesRef = {'na','na','"spike glycoprotein"','na','na', ...       % Some proteins are named differently
    '"membrane protein"','na','na','na','na','"nucleocapsid protein"','na'};

%% Needleman-Wunsch alignment of protein sequences.
for s = 1:length(seqs)
    proteinSequences2 = seqs(s).Translation.Sequences;                      % Sequences of the proteins present in virus s.
    proteinNames2 = seqs(s).Translation.Names;                              % Names of the corresponding proteins

    alignments = char();                                                    % Initialize an array with the alignments
    for p=1:length(proteinNames2)                                           % Loop through the proteins present
        protein = proteinNames2{p};                                         % Name of the current protein                                        
        % Find the corresponding protein in reference sequence:
        if (contains(protein,proteinNamesRef))==true | ...               	% If there is a corresponding protein in the reference           
                (contains(protein,proteinNamesRefLower))==true | ...        % (can also be in lower case letters), ...
                (contains(protein,proteinAlternativeNamesRef))==true        % (can also be one of the alternative names), ...
            m=1;
            while strcmp(protein, proteinNamesRef{m})==false & ...          % Find the index 'm' of the corresponding protein in the reference.
                    strcmp(protein, proteinNamesRefLower{m})==false & ...
                    strcmp(protein, proteinAlternativeNamesRef{m})==false
                m = m+1;
            end
            Seq1 = proteinSequencesRef{m};                                  % Reference sequence.
            Seq2 = proteinSequences2{p};                                    % Query sequence
            Seq2(find(Seq2=='J')) = 'I';                                    % Convert 'J' into 'I' (there is a bug in sequence 1324).
            [score, alignment, start] = nwalign(Seq1,Seq2);                 % Needleman-Wunsch alignment
            alignments{p} = alignment;
        else
            alignments{p} = [];
        end
    end
    seqs(s).Translation.Alignment = alignments;                             % Store the alignments of all proteins in seqs.Translation
    if mod(s,50)==0                                                         % Show progression.
        disp(s)
    end
end
%% Save results
save('Sequence_alignment.mat','names','seqs')