function [xMAP] = plotcorner (x, P, x0, xbnds, count, BurnIn, VarNames)
%
% [xMAP] = plotcorner (x, P, x0, xbnds, count, BurnIn, VarNames)
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


[~, Nvar] = size(x);

% find the maximum a posteriori model
[~, xMAPind] = max(P);
xMAP = x(xMAPind,:);

VarVary = diff(xbnds,[],2)>0;

if isempty(x0), x0 = x(1,:); end

figure;
figpos = get(0,'ScreenSize');
set(gcf,'Position',[figpos(1),100,0.8*figpos(3),0.85*figpos(4)]);
hAx = tight_subplot(Nvar,Nvar,0,0.04,0.04);

for mi = 1:Nvar
    for ni = 1:Nvar
        axes( hAx((mi-1)*Nvar + ni) );
        
        if mi==ni
            
            histogram(x(BurnIn:end,mi), 40, 'EdgeColor', 'none');   hold on;
            % bounds
            plot(xbnds(mi,1)*ones(1,2), ylim, 'r-');
            plot(xbnds(mi,2)*ones(1,2), ylim, 'r-');
            
            % starting model
            plot(x0(mi)*ones(1,2), ylim, 'k:');
            
            % maximum a posteriori model
            plot(xMAP(mi)*ones(1,2), ylim, 'r:');
            hold off;
            
            if VarVary(mi); xlim(xbnds(mi,:)); end
            xlabel(VarNames{mi});
            
            set(gca,'YTickLabel',[]);
            
        elseif mi<ni
            dscatter(x(BurnIn:end,ni), x(BurnIn:end,mi),'plottype','contour');
            set(gca,'XTickLabel',[]);
            if VarVary(ni); xlim(xbnds(ni,:)); end
            if VarVary(mi); ylim(xbnds(mi,:)); end
            if ni==Nvar
                set(gca,'YAxisLocation','right'); 
            else
                set(gca,'YTickLabel',[]);
            end
        else
            set(gca,'visible', 'off');
        end
    end
end

% sgtitle('Correlation Plot');
annotation('textbox','Position',[0.45,0.9,0.1,0.1],...
    'String','Correlation Plot (black line: starting model; red line: MAP model)',...
    'HorizontalAlignment','center','FontSize',20,'EdgeColor','none');
end