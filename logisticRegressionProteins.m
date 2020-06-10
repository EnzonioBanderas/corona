function [beta0, beta1] = logisticRegressionProteins()
% Function performs a logistic regression on the number of mutations for 
% protein. The mutation numbers are stored in 'mismatchBoolean'.
% beta0 and beta1 are arrays with the values of beta1 and beta2 for all
% proteins.

    % load the data
    data = load('Data/mismatchBoolean0602.mat');
    mismatchBoolean = data.mismatchBoolean;
    proteinNames = fields(mismatchBoolean);
    % remove NSP11:
    proteinNames(11) = [];
    % initialize beta arrays:
    beta0 = zeros(1,length(proteinNames));
    beta1 = zeros(1,length(proteinNames));
    % loop over proteins
    for i = 1:length(proteinNames)
        protein = proteinNames{i};
        y = mismatchBoolean.(protein)';
        X = linspace(0,length(y)-1,length(y))';
        % Define extra cost on first datapoint (x = 0)
        alpha = 5*max(find(y));
        if strcmp(protein, 'NSP3')
            alpha = 20000;
        end
        % Do logistic regression
        b = logisticRegression(X, y, alpha);
        beta0(i) = b(1);
        beta1(i) = b(2);
    end
end

%% Functions
function g = sigmoid(z)
% Calculate sigmoid of z.
    g = 1 ./ (1 + exp(-z));
end

function [J, grad] = computeCost(b, X, y, alpha)
% Calculate cost and gradient.
% Extra cost factor is put on the first term (x=0). 
    m = size(X, 1);
    X = [ones(m, 1) X];                                                     % add constant term
    
    h = sigmoid(X * b);
    
    J1 = -(alpha / m) * ( y(1) * log(h(1)) + (1 - y(1)) * log(1 - h(1)) );
    J2 = -(1 / m) * sum( y(2:end) .* log(h(2:end)) + (1 - y(2:end)) .* log(1 - h(2:end)));
    J = J1 + J2;
    
    grad = zeros(size(b,1), 1);
    
    for i = 1:size(grad)
        grad1 = (alpha / m) * (h(1) - y(1)) * X(1,i) ;
        grad2 = (1 / m) * sum( (h(2:end) - y(2:end))' * X(2:end,i) );
        grad(i) = grad1 + grad2;
    end
end

function b = logisticRegression(X, y, alpha)
    % Compute cost and gradient
    [m, n] = size(X);
    b0 = zeros((n+1),1);
    % Gradient descent
    options = optimset('GradObj', 'on', 'MaxIter', 400);
    b = fminunc(@(t)computeCost(t, X, y, alpha), b0, options);
end
