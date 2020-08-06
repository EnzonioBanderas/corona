function r = replicationRate(d, r0, sigma, beta)
% This function determines the fitness of a strain with a distance d from
% the reference sequence, which has fitness r0.
% The fitness is normal distributed with mean mu (predicted by logistic 
% regression) and standard deviation sigma.

    pNames = fields(pInfo);
    rArray = zeros(length(pNames),1);
    for i=1:length(pNames)
        if d(i) == 0
            replRate = 1;
        else
            b = beta(:, i);
            mu = predict(b, d(i));
            replRate = max(sigma(i)*randn + mu, 0);
        end
        rArray(i) = replRate;
    end
    r = r0 * prod(rArray);
    %r = rArray;
end