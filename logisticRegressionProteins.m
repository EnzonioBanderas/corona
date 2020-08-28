function [beta, sigma] = logisticRegressionProteins()
% Function performs a logistic regression on the number of mutations for 
% protein. The mutation numbers are stored in 'mismatchBoolean'.
% beta0 and beta1 are arrays with the values of beta1 and beta2 for all
% proteins.
% Sigma is an array with the standard error of the logistic regression for
% each protein

    % load the data
    data = load('Data/mismatchBoolean0821.mat');
    mismatchBoolean = data.mismatchBoolean;
    proteinNames = fields(mismatchBoolean);
    % initialize beta arrays:
    beta = zeros(2,length(proteinNames));
    sigma = zeros(1, length(proteinNames));
    % loop over proteins
    for i = 1:length(proteinNames)
        protein = proteinNames{i};
        y = mismatchBoolean.(protein)';
        N = length(y);
        X = transpose(0 : N-1);
 
        % Do logistic regression
        b = logisticRegression(X, y);
        beta(:, i) = b;
        
        % Find standard error
        prediction = predict(b, X);
        sigma(i) = sqrt (sum((prediction - y) .^2) / (N - 2) );
    end
end

%% Functions
function g = sigmoid(z)
% Calculate sigmoid of z.
    g = 1 ./ (1 + exp(-z));
end

function [J, grad] = computeCost(b, X, y)
% Calculate cost and gradient.
% Extra cost factor is put on the first term (x=0). 
    m = size(X, 1);
    % add constant term
    X = [ones(m, 1) X];                                                    
    
    h = sigmoid(X * b);
    
    J = -(1 / m) * sum( y .* log(h) + (1 - y) .* log(1 - h));
    
    grad = zeros(size(b,1), 1);
    
    for i = 1:size(grad)
        grad(i) = (1 / m) * sum( (h - y)' * X(:,i) );
    end
end

function b = logisticRegression(X, y)
    % Compute cost and gradient
    [m, n] = size(X);
    b0 = zeros((n+1),1);
    % Gradient descent
    options = optimset('GradObj', 'on', 'MaxIter', 400, 'Display', 'off');
    b = fminunc(@(t)computeCost(t, X, y), b0, options);
end

