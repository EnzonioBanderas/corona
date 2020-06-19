clear all
close all

maxd = 9;
dArray = linspace(0,maxd,maxd+1);
L = 30000;
%%
sigma = zeros(1,length(dArray));
for i = 1:length(dArray)
    d = dArray(i);
    sigmaArray = zeros(1,d);
    for j = 1:d
        sigmaArray(j) = L - j;
    end
    sigma(i) = (3^d) * prod(sigmaArray);
end

semilogy(dArray, sigma)

function n_factorial = stirling(n)
    n_factorial = sqrt(2*pi.*n) .* (n ./ exp(1)).^n;
end

% clear all
% close all
% 
% [gRefSeq, pRefSeq, pNames, pInfo] = getRefSeq();
% [beta0, beta1] = logisticRegressionProteins();
% for i = 1:length(pNames)
%     protein = pNames{i};
%     pInfo.(protein).betas(1) = beta0(i);
%     pInfo.(protein).betas(2) = beta1(i);
% end
% %%
% r0 = 1;
% Nmax = 40;                  % test 40 distances
% Nproteins = length(pNames);
% sigmas = linspace(0.01, 0.1, 4);
% Nrep = 200;
% figure()
% for s = 1:length(sigmas)
%     sigma = sigmas(s);
%     r = zeros(Nrep, Nmax+1, Nproteins);
%     for n = 1:Nrep
%         
%         d = repmat(linspace(0,Nmax,Nmax+1),Nproteins,1);
%         
% 
%         for i=1:Nmax+1
%             r(n,i,:) = replicationRate(d(:,i), r0, sigma, pInfo);
%         end
%     end
% 
%     subplot(2,2,s)
%     for i = 1:1
%         p1 = plot(d(i,:), r(:,:,i),'.r');
%         hold on
%         protein = pNames{i};
%         b = pInfo.(protein).betas;
%         mu = predict(b', d(i,:)');
%         p2 = plot(d(i,:), mu, '-b','LineWidth',2);
%         ylim([0,1.4])
%     end
%     hold off
%     set(gca, 'LineWidth',2,'FontSize',20)
%     xlabel('Number of mutations')
%     ylabel('Relative reproduction rate')
%     %legend('show')
%     legend([p1(1), p2],['Simulated reproduction rate (N=', num2str(Nrep)], 'Mean reproduction rate (logistic curve)')
%     title(['Sigma = ',num2str(sigma)])
% end
% suptitle('Simulated reproduction rate for different sigmas')