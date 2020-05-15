%{
Find country and time of modification.
%}

%% Load fasta file
file = 'sequences.fasta.txt';
data = fastaread(file);   

%% Time and location analysis
j=0;
for i = 1:length(data)
    L = length(data(i).Sequence);                                           % length of sequence
    N = count(data(i).Sequence, 'N');                                       % number of unresolved nucleotides ('N' in fasta file)
    if L-N > 28000 & L-N < 30000                                            % check if the number of unresolved nucleotides (L-N) is between 28000 and 30000.
        j = j+1;
        header = split(data(i).Header, '|');
        headers{j} = header;
        presentLocs{j} = header{4};
    end
end
uniqueLoc = unique(presentLocs);
wrongLoc = uniqueLoc(1:11);

j=0;
for i = 1:length(data)
    L = length(data(i).Sequence);                                           % length of sequence
    N = count(data(i).Sequence, 'N');                                       % number of unresolved nucleotides ('N' in fasta file)
    if L-N > 28000 & L-N < 30000                                            % check if the number of unresolved nucleotides (L-N) is between 28000 and 30000.
        j = j+1;
        cutData(j) = data(i);
        header = split(data(i).Header, '|');
        if ismember(header{4},wrongLoc)
            N = length(header);
            locs{j} = header{N-3};
        else
            locs{j} = header{4};
        end
        header2 = split(data(i).Header, ' |');
        name = header2(1);
        if strcmp(name,seqs(j).LocusName)==1
            seqs(j).('Country') = locs{j};
        else
            disp(j)
        end
    end
end

[~,~,ic] = unique(locs);                                       % Unique Values By Row, Retaining Original Order
h = accumarray(ic, 1);                                                      % Count Occurrences
counts = h(ic);                                                             % Map Occurrences To ‘ic’ Values
uniqueLocs = cell(length(locs),1);
for i=1:length(locs)
    uniqueLocs{i} = join([locs{i},' ',num2str(counts(i))]);
end
final_uniqueLocs = unique(uniqueLocs);

%% Location analysis
Europe = {'Europe','Czech Republic','Finland','France','Germany','Greece','Italy',...
    'Netherlands','Spain','Sweden','Turkey'};
Asia = {'Asia','China','Hong Kong','India','Iran','Israel','Japan','Kazakhstan',...
    'Nepal','Pakistan','South Korea','Sri Lanka','Taiwan'};
SouthAmerica = {'SouthAmerica','Brazil','Colombia','Peru','Puerto Rico'};
NorthAmerica = {'NorthAmerica','USA'};
Oceania = {'Oceania','Australia','Malaysia','Thailand','Viet Nam'};
Africa = {'Africa','South Africa'};

continents = {Africa,Asia,Europe,SouthAmerica,NorthAmerica,Oceania};
% Find continent
for i = 1:length(seqs)
    for c = 1:length(continents)
        if ismember(seqs(i).Country, continents{c})
            seqs(i).Continent = continents{c}{1};
        end
    end
end

%% Time analysis
for i=1:length(seqs)
    seqs(i).LocusModificationDate = datetime(seqs(i).LocusModificationDate);
end
%%
for i=1:length(seqs)
    dates(i) = seqs(i).LocusModificationDate;
end
%%
periods = {'FEB','MARCH','APR1','APR2','MAY'};
for i=1:length(seqs)
    date = seqs(i).LocusModificationDate;
    if date < datetime('01-MAR-2020')
        seqs(i).Period = 'FEB';
    elseif date < datetime('01-APR-2020') & date >= datetime('01-MAR-2020')
        seqs(i).Period = 'MARCH';
    elseif date < datetime('16-APR-2020') & date >= datetime('01-APR-2020')
        seqs(i).Period = 'APR1';
    elseif date < datetime('01-MAY-2020') & date >= datetime('16-APR-2020')
        seqs(i).Period = 'APR2';
    elseif date >= datetime('01-MAY-2020')
        seqs(i).Period = 'MAY';
    end
end