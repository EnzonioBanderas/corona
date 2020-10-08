% latin hypercube intitialization

% lh initial test, 1e3 points (probably 1e1 iterations per point, so 1e4 in total)
params_lhs = lhsdesign(1e3, 7, 'Iterations', 1e4); % [a b r0 mu]
save('Data/params_lhs_1e3_7_Iter1e4.mat', 'params_lhs')

% lh initial test, 1e4 points (probably 1e1 iterations per point, so 1e5 in total)
params_lhs = lhsdesign(1e4, 7, 'Iterations', 1e3); % [a b r0 mu]
save('Data/params_lhs_1e4_7_Iter1e3.mat', 'params_lhs')

% lh initial test, 1e5 points (probably 1e1 iterations per point, so 1e6 in total)
params_lhs = lhsdesign(1e5, 7, 'Iterations', 1e2); % [a b r0 mu]
save('Data/params_lhs_1e5_7_Iter1e2.mat', 'params_lhs')