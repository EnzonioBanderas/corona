clearvars

% Define
possible_bases = ['A', 'C', 'G', 'T'];
nPossibleBases = length(possible_bases);

% Initial conditions
sequence_initial = {'AAAAAAAAAA', 'GGGGGGGGGG'}; % sequence(s) of initial virus(es) (these are deemed "viable")
sequence = sequence_initial; % represents sequences which are currently present
nSequenceInitial = length(sequence); % number of initial sequences
sequence_number = [1, 2]; % initial virus count(s)
sequence_growth = [0.2, 0.3]; % probability of growth for each sequence
sequence_mu_v = [0.01, 0.02]; % for every replication, every base has a chance to change from its original with this probability
sequence_recognition = false(1, nSequenceInitial); % initial virus is not recognized
nBase = length(sequence{1}); % number of bases is constant

% Parameters - probabilities and number of time ticks
growth_mutation_range = 0.3; % range in which growth can change in a single mutation
mu_v_mutation_range = 0.3; % range in which mu_v can change in a single mutation
immune_recognition = 0.001; % each time tick, this is the probability that a virus will be recognized for every time tick
immune_destruction = 0.75; % if the virus is recognized, this is the probability that a virus will be eliminated for every time tick
nTime = 40; % number of time ticks
viability = ones(1, nBase); % viability of viral sequence as a function of number of mutations
viability(5:end)=0; % sequences with 4 or more mutations from sequences which are known to be viable initially are destroyed

% Initialization of variables which record viral info over time
t=0:nTime;
recorded_sequence=sequence;
recorded_sequence_time=cell(1, nSequenceInitial);
recorded_sequence_growth=cell(1, nSequenceInitial);
recorded_sequence_mutation=cell(1, nSequenceInitial);
recorded_sequence_number=cell(1, nSequenceInitial);
recorded_sequence_recognition_time=zeros(1, nSequenceInitial);
recorded_sequence_recognition=false(1, nSequenceInitial);
for iSequence=1:nSequenceInitial
    recorded_sequence_time{iSequence}(1)=0;
    recorded_sequence_number{iSequence}(1)=sequence_number(iSequence);
    recorded_sequence_growth{iSequence}(1)=sequence_growth(iSequence);
    recorded_sequence_mutation{iSequence}(1)=sequence_mu_v(iSequence);
end
average_mutation = zeros(1, nTime+1);
average_mutation(1) = sum(sequence_number.*sequence_mu_v)/sum(sequence_number);
average_growth = zeros(1, nTime+1);
average_growth(1) = sum(sequence_number.*sequence_growth)/sum(sequence_number);
sequence_viable_initial = true(1, nSequenceInitial); % viability when compared to initial sequences

for iTime=1:nTime
    
    % Go through each different virus and evaluate whether it grows or mutates 
    nSequence = length(sequence); % number of different sequences
    for iSequence=1:nSequence % current sequence is sequence{iSequence}
        
        initialNumber = sequence_number(iSequence); % number of viruses
        for iSequenceNumber=1:initialNumber % current virus
            
            % Virus can grow depending on its growth probability
            if rand(1)<sequence_growth(iSequence)

                % When the virus replicates, it might mutate. The mutation
                % probability is applied on each base.
                new_sequence=sequence{iSequence}; % initialize new sequence as old sequence
                for iBase=1:nBase
                    
                    % if a mutation happens, change the base to another one
                    if rand(1)<sequence_mu_v(iSequence)
                        current_base = sequence{iSequence}(iBase);
                        current_possible_bases = possible_bases(possible_bases~=current_base);
                        new_base = current_possible_bases(ceil(rand(1)*(nPossibleBases-1)));
                        new_sequence(iBase) = new_base;
                    end
                    
                end

            % if the new sequence generated is part of old sequences, 
            % simply add 1 to the respective sequence virus number
            new_sequence_comparison = strcmp(sequence, new_sequence);
            if any(new_sequence_comparison)
                new_sequence_index = find(new_sequence_comparison);
                sequence_number(new_sequence_index)=sequence_number(new_sequence_index)+1;
                
            % else it must be a new sequence
            else
                
                % Check whether new sequence is viable compared to all
                % initial sequences
                for iInitialSequence=1:nSequenceInitial
                    
                    % Number of mutations (distance from initial sequences)
                    nMutation = sum(new_sequence~=sequence_initial{iInitialSequence});
                    
                    % Check viability as a function of number of mutations
                    % and assign boolean for each initial sequence
                    if nMutation>0
                        sequence_viable_initial(iInitialSequence) = viability(nMutation); 
                    end
                    
                end
                
                % If the sequence is viable compared to any of the
                % initial sequences, it is viable.
                % The new sequence can be added to the sequence cell array.
                if any(sequence_viable_initial)
                    
                    % Add new sequence
                    nNewSequence=length(sequence)+1;
                    sequence{nNewSequence}=new_sequence;
                    sequence_number(nNewSequence)=1;
                    sequence_recognition(nNewSequence)=false;                

                    % Add new sequence growth probability by randomly changing
                    % the growth probability of the unmutated sequence in a
                    % defined range.
                    sequence_growth(nNewSequence)=sequence_growth(iSequence)+(rand(1)-0.5)*growth_mutation_range;
                    if sequence_growth(nNewSequence)<0
                        sequence_growth(nNewSequence)=0;
                    elseif sequence_growth(nNewSequence)>1
                        sequence_growth(nNewSequence)=1;
                    end

                    % Add new sequence mutation probability by randomly changing
                    % the mutation probability of the unmutated sequence in a
                    % defined range.
                    sequence_mu_v(nNewSequence)=sequence_mu_v(iSequence)+(rand(1)-0.5)*mu_v_mutation_range;
                    if sequence_mu_v(nNewSequence)<0
                        sequence_mu_v(nNewSequence)=0;
                    elseif sequence_mu_v(nNewSequence)>1
                        sequence_mu_v(nNewSequence)=1;
                    end
                
                end
                
            end
            
            end
        end
    end
    
    % Evaluate for each unrecognized sequence whether it is recognized
    sequence_unrecognized = ~sequence_recognition;
    for iUnrecognized=find(sequence_unrecognized)
        
        % if a sequence was already recorded, assign its recorded recognition value
        sequence_in_recorded_sequence = strcmp(recorded_sequence, sequence{iUnrecognized});
        if recorded_sequence_recognition(sequence_in_recorded_sequence)
                sequence_recognition(iUnrecognized) = true; 
                continue % recognized, so continue to next unrecognized
        end
        
        for iSequenceNumber=1:sequence_number(iUnrecognized)
            % add here that if a sequence was already recognized previously
            % it should be immediately set to recognized again
            if rand(1)<immune_recognition
                sequence_recognition(iUnrecognized) = true; 
                recorded_sequence_recognition_time(strcmp(recorded_sequence, sequence(iUnrecognized))) = t(iTime);
                break % recognized, so continue to next unrecognized
            end
        end
    end
    
    % For recognized viruses, evaluate whether they are destroyed
    for iRecognized=find(sequence_recognition)
        initialNumber=sequence_number(iRecognized);
        for iSequenceNumber=1:initialNumber
            if rand(1)<immune_destruction
                sequence_number(iRecognized)=sequence_number(iRecognized)-1;
            end
        end
    end
    
    % Add info over time for viruses
    recorded_sequence=unique([recorded_sequence, sequence], 'stable');
    recorded_sequence_time=[recorded_sequence_time, cell(1, length(recorded_sequence)-length(recorded_sequence_time))];
    recorded_sequence_mutation=[recorded_sequence_mutation, cell(1, length(recorded_sequence)-length(recorded_sequence_mutation))];
    recorded_sequence_growth=[recorded_sequence_growth, cell(1, length(recorded_sequence)-length(recorded_sequence_growth))];
    recorded_sequence_number=[recorded_sequence_number, cell(1, length(recorded_sequence)-length(recorded_sequence_number))];
    recorded_sequence_recognition=[recorded_sequence_recognition, false(1, length(recorded_sequence)-length(recorded_sequence_number))];
    for iRecordedSequence=1:length(recorded_sequence)
        
        recorded_sequence_in_sequence = strcmp(sequence, recorded_sequence{iRecordedSequence});
        
        if any(recorded_sequence_in_sequence)
            recorded_sequence_time{iRecordedSequence}=[recorded_sequence_time{iRecordedSequence}, iTime];
            recorded_sequence_mutation{iRecordedSequence}=[recorded_sequence_mutation{iRecordedSequence}, sequence_mu_v(recorded_sequence_in_sequence)];
            recorded_sequence_growth{iRecordedSequence}=[recorded_sequence_growth{iRecordedSequence}, sequence_growth(recorded_sequence_in_sequence)];
            recorded_sequence_number{iRecordedSequence}=[recorded_sequence_number{iRecordedSequence}, sequence_number(recorded_sequence_in_sequence)];
            recorded_sequence_recognition(iRecordedSequence)=sequence_recognition(recorded_sequence_in_sequence); % also updating for old sequences...
        end
        
    end
    
    % Remove any sequences that have a sequence number of 0
    sequence_remaining = sequence_number~=0;
    sequence=sequence(sequence_remaining);
    sequence_number=sequence_number(sequence_remaining);
    sequence_recognition=sequence_recognition(sequence_remaining);
    sequence_growth=sequence_growth(sequence_remaining);
    sequence_mu_v=sequence_mu_v(sequence_remaining);
    
    % Update average mutation and growth
    average_mutation(iTime+1) = sum(sequence_number.*sequence_mu_v)/sum(sequence_number);
    average_growth(iTime+1) = sum(sequence_number.*sequence_growth)/sum(sequence_number);

end

% Display all generated sequences
disp(recorded_sequence)

% Plot time courses of virus numbers for each sequence
figure
nSequenceRecorded = length(recorded_sequence_number);
total_virus_number_end=0;
virus_number_end = zeros(1, nSequenceRecorded);
for iSequence=1:nSequenceRecorded
    plot(recorded_sequence_time{iSequence}, recorded_sequence_number{iSequence})
    hold on
    
    virus_number_end(iSequence) = recorded_sequence_number{iSequence}(end);
    total_virus_number_end=total_virus_number_end+recorded_sequence_number{iSequence}(end);
end
title('Virus count time series')
xlabel('Time')
ylabel('Virus count')

% Display total virus count
disp(total_virus_number_end)

% Plot sequences in 2D space by calculating distances between sequences and
% using dimensional reduction (MDS). Size and color of sequences varies 
% with virus number.
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
tic
recorded_sequence_distances = seqpdist(recorded_sequence,...
                                        'Alphabet', 'NT',...
                                        'Method', 'Jukes-Cantor',...
                                        'SQUAREFORM', true);
toc
tic
coord_MDS = cmdscale(recorded_sequence_distances, 2);
toc
scatter(coord_MDS(:, 1), coord_MDS(:, 2), 10*(virus_number_end+1), virus_number_end, 'filled')
title('Sequence scatter plot, dimensional reduction with MDS on Jukes-Cantor distance matrix')
xlabel('MDS 1')
ylabel('MDS 2')

% Plot average growth and mutation probabilities over time
figure
plot(t, average_growth)
hold on
plot(t, average_mutation)
title('Virus average mutation and growth probabilities')
xlabel('Time')
ylabel('Probability')
legend({'growth', 'mutation'})