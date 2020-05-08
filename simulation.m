clearvars

sequence = {'AAAAAAAAAA', 'GGGGGGGGGG'}; % sequence(s) of initial virus(es)
sequence_number = [1, 2]; % initial virus count(s)
nSequence = length(sequence); % number of viruses
sequence_recognition = false(1, nSequence); % initial virus is not recognized
possible_bases = ['A', 'C', 'G', 'T'];
nPossibleBases=length(possible_bases);
nBase = length(sequence{1}); % no change in length
growth_probability = 0.2; % every 5 time ticks I want a replication
mu_v = 0.1; % for every replication, this is the chance of a mutation to another base than original
nTime = 50; % number of time ticks
immune_recognition = 0.01; % each time tick, this is the probability a virus will be recognized
immune_destruction = 0.3; % if the virus is recognized, this is the chance that a virus will be eliminated

% time info initialization
t=0:nTime;
recorded_sequence=sequence;
recorded_sequence_time=cell(1, nSequence);
recorded_sequence_number=cell(1, nSequence);
recorded_sequence_recognition_time=zeros(1, nSequence);
recorded_sequence_recognition=false(1, nSequence);
for iSequence=1:nSequence
    recorded_sequence_time{iSequence}(1)=0;
    recorded_sequence_number{iSequence}(1)=sequence_number(iSequence);
end

for iTime=1:nTime
    
    % Evaluate for each virus whether it grows
    nSequence = length(sequence); % number of viruses
    for iSequence=1:nSequence
        initialNumber=sequence_number(iSequence);
        for iSequenceNumber=1:initialNumber
            if rand(1)<growth_probability

                % go through each base and evaluate whether it mutates
                new_sequence=sequence{iSequence};
                for iBase=1:nBase
                    % if a mutation happens, change the base to another one
                    if rand(1)<mu_v
                        current_base = sequence{iSequence}(iBase);
                        current_possible_bases = possible_bases(possible_bases~=current_base);
                        new_base = current_possible_bases(ceil(rand(1)*(nPossibleBases-1)));
                        new_sequence(iBase) = new_base;
                    end
                end

            % if the new sequence generated is part of old
            % sequences, add a number to that old sequence
            new_sequence_comparison = strcmp(sequence, new_sequence);
            if any(new_sequence_comparison)
                new_sequence_index = find(new_sequence_comparison);
                sequence_number(new_sequence_index)=sequence_number(new_sequence_index)+1;
                
            % else it must be new and you should grow the sequence array
            else
                nNewSequence=length(sequence)+1;
                sequence{nNewSequence}=new_sequence;
                sequence_number(nNewSequence)=1;
                sequence_recognition(nNewSequence)=false;
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
    recorded_sequence_number=[recorded_sequence_number, cell(1, length(recorded_sequence)-length(recorded_sequence_number))];
    recorded_sequence_recognition=[recorded_sequence_recognition, false(1, length(recorded_sequence)-length(recorded_sequence_number))];
    % length of recorded_sequence_time and recorded_sequence_number should
    % be the same as length of recorded_sequence, this can be done by
    % initialization
    for iRecordedSequence=1:length(recorded_sequence)
        
        recorded_sequence_in_sequence = strcmp(sequence, recorded_sequence{iRecordedSequence});
        
        if any(recorded_sequence_in_sequence)
            recorded_sequence_time{iRecordedSequence}=[recorded_sequence_time{iRecordedSequence}, iTime];
            recorded_sequence_number{iRecordedSequence}=[recorded_sequence_number{iRecordedSequence}, sequence_number(recorded_sequence_in_sequence)];
            recorded_sequence_recognition(iRecordedSequence)=sequence_recognition(recorded_sequence_in_sequence); % also updating for old sequences...
        end
        
    end
    
    % Remove any sequences that have a sequence number of 0
    sequence_remaining = sequence_number~=0;
    sequence=sequence(sequence_remaining);
    sequence_number=sequence_number(sequence_remaining);
    sequence_recognition=sequence_recognition(sequence_remaining);

end

recorded_sequence

figure
total_virus_number_end=0;
for iSequence=1:length(recorded_sequence_number)
    plot(recorded_sequence_time{iSequence}, recorded_sequence_number{iSequence})
    hold on
    
    total_virus_number_end=total_virus_number_end+recorded_sequence_number{iSequence}(end);
end
xlabel('Time')
ylabel('Virus count')

total_virus_number_end