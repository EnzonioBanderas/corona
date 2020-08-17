function genomeRefCell = genome2cell(genomeData)
% This function takes as input the entire SARS-CoV2 reference genome. It
% stores the sequences of 24 genes encoding for the 24 proteins in a cell array.
    genomeRefCell{1} = genomeData(266-3:805+3);
    genomeRefCell{2} = genomeData(806-3:2719+3);
    genomeRefCell{3} = genomeData(2720-3:8554+3);
    genomeRefCell{4} = genomeData(8555-3:10054+3);
    genomeRefCell{5} = genomeData(10055-3:10972+3);
    genomeRefCell{6} = genomeData(10973-3:11842+3);
    genomeRefCell{7} = genomeData(11843-3:12091+3);
    genomeRefCell{8} = genomeData(12092-3:12685+3);
	genomeRefCell{9} = genomeData(12686-3:13024+3);
    genomeRefCell{10} = genomeData(13025-3:13441+3);
    genomeRefCell{11} = [genomeData(13442-3:13467+3), genomeData(13467-3:16236+3)];
    genomeRefCell{12} = genomeData(16237-3:18039+3);
    genomeRefCell{13} = genomeData(18040-3:19620+3);
    genomeRefCell{14} = genomeData(19621-3:20658+3);
    genomeRefCell{15} = genomeData(20659-3:21552+3);
    genomeRefCell{16} = genomeData(21563-3:25381+3); % Spike
	genomeRefCell{17} = genomeData(25393-3:26217+3); % NS3
    genomeRefCell{18} = genomeData(26245-3:26469+3); % E
    genomeRefCell{19} = genomeData(26523-3:27188+3); % M
    genomeRefCell{20} = genomeData(27202-3:27384+3); % NS6
    genomeRefCell{21} = genomeData(27394-3:27756+3); % NS7a
    genomeRefCell{22} = genomeData(27756-3:27884+3); % NS7b
    genomeRefCell{23} = genomeData(27894-3:28256+3); % NS8
    genomeRefCell{24} = genomeData(28274-3:29530+3); % N
end