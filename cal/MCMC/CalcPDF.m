function [m, prob] = CalcPDF (mbnds, models, Nbins)
% [m, prob] = CalcPDF (mbnds, models, Nbins)
% 
% converts ensemble of models into a probability density function with
% points specified at m. 

if isempty(mbnds), mbnds = [min(models), max(models)]; end

[Niter, Nvar] = size(models);
m    = zeros(Nbins, Nvar); % bin centers
prob = zeros(Nbins, Nvar); % probabilities

for mi = 1:Nvar
    
    if mbnds(mi,1)==mbnds(mi,2), continue; end
    
    dm = (mbnds(mi,2)-mbnds(mi,1))/(Nbins-1);
    edges = (mbnds(mi,1)-0.5*dm):dm:(mbnds(mi,2)+0.5*dm);
    m(:,mi) = mbnds(mi,1):dm:mbnds(mi,2);
    
    counts = histcounts(models(:,mi),edges);
    prob(:,mi) = counts/Niter/dm;
    
end


end