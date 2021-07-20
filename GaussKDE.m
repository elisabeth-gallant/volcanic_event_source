function [GMM, kernelBW] = GaussKDE(vents)
    Nv = size(vents,2);
    C = cov(vents');
    [~,~,V] = svd(C);
    basis = V*(vents-mean(vents,2));
    [~, ~, basis] = ksdensity(basis');
    kernelBW = V*diag(basis.^2)*V';
    GMM = gmdistribution(vents',repmat(kernelBW,[1,1,Nv]));
end