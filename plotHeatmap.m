clear all

% filePath = dir('data_iteration/data_iteration*');
% nFiles = length(filePath)
% data_cell = cell(1, nFiles)
% for iFile = 1:nFile
%   load(file(iFile))
%   data_cell{iFile} = data
%   params_cell{iFile} = params
% end
%
% get responses and input by concatenation, first between structures and
% then between structure array fields

% data_concat = [data_cell{:}]

% [data_concat(:).c]
% vertcat(data_concat(:).c)

titer_fn = dir('titer*/titer*');
load([titer_fn(1).folder, '\', titer_fn(1).name])
nTiter_fn = length(titer_fn);
la(1:nTiter_fn)=titer;
num = zeros(1, nTiter_fn);
Mu = num;
a = num;
U_sum = num;
statY = nan(1, nTiter_fn);
statD = statY;
statR = statY;
for iIter = 1:nTiter_fn
    load([titer_fn(iIter).folder, '\', titer_fn(iIter).name])
%     titer_cat(iIter) = titer;

    % add data in titer which is always present to separate arrays
    Mu(iIter) = titer.Mu;
    a(iIter) = titer.a;
    U_sum(iIter) = titer.U_sum;

    % if statY exists in titer then statR and statD also exist and should
    % be added to separate arrays
    if isfield(titer, 'statY')
        statY(iIter) = titer.statY;
        statD(iIter) = titer.statD;
        statR(iIter) = titer.statR;
    end
    
    % Get seed so that you can possible later filter out duplicate seeds
    name_split = split(titer_fn(iIter).name, '_');
    name_split = split(name_split{2}, '.');
    num(iIter) = str2double(name_split{1});
end



% X = [titer_cat.Mu];
% Y = [titer_cat.a];
% V = [titer_cat.U_sum];




N = 100;
% xq = logspace(log10(min(X)), log10(max(X)), N);
% xq = logspace(-6, -3, N+1);
Muq = logspace(-6, -3, N);
Muq_i = 1:N;
% yq = linspace(min(Y), max(Y), N);
% yq = linspace(1e-5, 3e-5, N+1);
aq = linspace(1e-5, 3e-5, N);
aq_i = Muq_i;



%% Plot U_sum, statD, statY and statR
data = {U_sum, statD, statY, statR};
data_str = {'Number of uninfected cells' , 'statD' , 'statY' , 'statR'};

for iData = 1:length(data)
    
    [MUq, Aq] = meshgrid(Muq, aq);
    MUq = reshape(MUq, [1, N^2]);
    Aq = reshape(Aq, [1, N^2]);

    % Vq = gridfit(X,Y,V,xq,yq);
    isnotnan = ~isnan(data{iData});
    data_plot = data{iData}(isnotnan);
    F = scatteredInterpolant(Mu(isnotnan)', a(isnotnan)', data_plot', 'linear');
    data_plot = F(MUq', Aq');

    MUq = reshape(MUq', [N, N]);
    Aq = reshape(Aq', [N, N]);
    data_plot = reshape(data_plot, [N, N]);

    figure, imagesc(Muq, aq, data_plot)
    set(gca, 'XScale', 'log')
    set(gca, 'YDir', 'normal')
%     set(gca, 'ColorScale', 'linear')
    cb = colorbar;
    cb.Label.String = data_str{iData};
    xlabel('Mutation rate')
    ylabel('Infection rate')

end



% % Average values per bin, possibly do linear interpolation to fill in nan
% % values
% [NN,Xedges,Yedges,binX,binY] = histcounts2(X, Y, xq, yq);
% [Uniq, G, GG] = unique([binX; binY]', 'rows');
% V_mean = splitapply(@mean, V',  GG);
% Vq = nan(N);
% Vq(sub2ind([N, N], Uniq(:,2), Uniq(:,1))) = V_mean;