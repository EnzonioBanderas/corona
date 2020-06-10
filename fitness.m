function r = fitness(d, r0, sigma, pInfo)
% This function determines the fitness of a strain with a distance d from
% the reference sequence, which has fitness r0.
% The fitness is normal distributed with mean mu (predicted by logistic 
% regression) and standard deviation sigma.

    pNames = fields(pInfo);
    rArray = zeros(length(pNames),1);
    for i=1:length(pNames)
        protein = pNames{i};
        b = pInfo.(protein).betas;
        mu = predict(b', d(i));
        rArray(i) = max(sigma*randn + mu, 0);
    end

    %r = rArray;
    r = r0 * prod(rArray);
end