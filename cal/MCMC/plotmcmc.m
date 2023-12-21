function [xMAP] = plotmcmc (x, P, x0, xbnds, anneal, VarNames)
%
% [xMAP] = PlotMCMCAnalytics (x, P, x0, xbnds, count, BurnIn, VarNames)
% plots the outputs of mcmc inversion
%
% INPUTS
% x         matrix of model parameters from mcmc [Niter x Nvar]
% P         vector of posterior probabilities [Niter x 1]
% x0        starting model [1 x Nvar, can be empty]
% xbnds     bnds of parameters [Nvar x 2]
% count     number of accepted models [scalar]
% BurnIn    burn in of markov chain [scalar]
% VarNames  name of variables
%
% OUTPUT
% xMAP      maximum a posteriori model [1 x Nvar]
%

ivar = find(diff(xbnds,1,2)>0);
xvar = x(:,ivar);
[Niter, Nvar] = size(xvar);

Nrow = floor(sqrt(Nvar));
Ncol = ceil(Nvar/Nrow);

% find the maximum a posteriori model
[~, xMAPind] = max(P(anneal.burnin:end));
xMAP = x(xMAPind+anneal.burnin-1,:);

if isempty(x0), x0 = x(1,:); end

anneal.ilvl = [anneal.burnin,round(anneal.burnin+(1:1:anneal.levels)*(Niter-anneal.burnin-anneal.refine)/anneal.levels)];

figure(101);

% plot distributions of the model parameters
for mi = 1:Nvar
    subplot(Nrow,Ncol,mi);
    
    % plot the distribution
    histogram(x(1:anneal.burnin,ivar(mi)), min(Niter/10,100), 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none','Normalization','count'); hold on;
    for i=1:anneal.levels
        histogram(x(anneal.ilvl(i):anneal.ilvl(i+1),ivar(mi)), min(Niter/10,100), 'FaceColor',(1-(i-1)/(anneal.levels+1))*[0.15 0.05 0.9], 'EdgeColor', 'none','Normalization','count');
    end

    % bounds
    plot(xbnds(ivar(mi),1)*ones(1,2), ylim, 'r-');
    plot(xbnds(ivar(mi),2)*ones(1,2), ylim, 'r-');
    
    % starting model
    plot(x0(ivar(mi))*ones(1,2), ylim, 'k:');
    
    % maximum a posteriori model
    plot(xMAP(ivar(mi))*ones(1,2), ylim, ':','Color',[0.6 0.0 0.1]); 
    axis tight; box on; hold off;
    
    title(VarNames{ivar(mi)});
end
sgtitle('Posterior PDFs (black: starting model; red: best-fit)');


end