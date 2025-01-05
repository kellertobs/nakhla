function [x_keep,P_keep,count,xbest] = mcmc(dhatFunc,PriorFunc,LikeFunc,ConstrFunc,ConstrInd,x0,xbnds,anneal,Niter)
% 
% [x_keep,P_keep,count] = mcmc(dhatFunc,PriorFunc,LikeFunc,x0,xstep,xbnds,Niter)
% 
% performs Markov Chain Monte Carlo sampling using random walk
% and Metropolis Hastings. Allows any sort of distribution for prior and
% likelihood
% Ying-Qi Wong, 3 Sep 2019
%
% INPUTS:
% dhatFunc  = function that predicts data from model params
% priorFunc = function that evaluates log prior probabilities 
% likeFunc  = function that calculates log data likelihood
% x0        = initial estimate of parameter vector (Nvar x 1)
% xstep     = step size in all parameter directions (Nvar x 1)
% xbnds     = Nx2 matrix of lower and upper bounds (can be empty)
% Niter     = number of iterations
% 
% NB: priorFunc and likeFunc should only take one input each: the model
% parameter vector and the data vector respectively.
%  
% OUTPUTS:
% x_keep = array of samples (Nvar x Niter)
% P_keep = posterior distribution
% count  = number of accepted.  Acceptance ratio is count/Niter


%Analyze inputs
Nvar  = length(x0);          % find number of model parameters
x1    = x0;                  % set initial guess as first candidate
xbest = x1;
if isempty(xbnds), xbnds = [-Inf, Inf]; end

% check functions work
dhatFunc  = fcnchk(dhatFunc);         
PriorFunc = fcnchk(PriorFunc);
LikeFunc  = fcnchk(LikeFunc);

% evaluate first model
dhat  = dhatFunc(x1);    
Px1   = PriorFunc(x1); 
Ld_x1 = LikeFunc(dhat,x1); 
Px1_d = Px1 + Ld_x1; 
Pbest = Px1_d*10;
ibest = anneal.burnin;

%Initialize the vectors x_keep and L_keep to store the accepted models
count  = 0;
x_keep = zeros(Nvar,Niter);
P_keep = zeros(Niter,1);
Nprint = min(100,floor(Niter/100));
anneal.hist = zeros(Niter,1);
anneal.ilvl = round(anneal.burnin+(1:1:anneal.levels)*(Niter-anneal.burnin-anneal.refine)/anneal.levels);
anneal.temp = 1;
anneal.lvl  = 1;
AcceptRatio = nan;

fprintf(1,'\n\nRunning MCMC...\n\n');
progressbar('_start');


%Begin loop to perform MCMC
for i=1:Niter
    
    % print the progress
    if mod(i,Nprint)==0, progressbar(i/Niter*100); end
    
    %Get step size with annealing
    if any(i==anneal.ilvl)
        anneal.temp = anneal.temp/2;
    elseif i>Niter-anneal.refine
        anneal.temp = 1/2^anneal.levels;
    end

    xstep  = anneal.initstep .* anneal.temp;

    flag = 0;
    
    while flag == 0
        %Random walk chain to find the proposed model x2
%         x2 = anneal.temp.*x1 + (1-anneal.temp).*xbest + xstep.*(2*rand(Nvar,1)-1);
        if i>Niter-anneal.refine; x1 = xbest; end
        x2 = x1 + xstep.*(2*rand(Nvar,1)-1);

        %Apply specified constraints (e.g., sum constraints, etc.)
        x2(ConstrInd) = ConstrFunc(x2(ConstrInd));

        %Check that the proposed model falls within the bounds.  If it falls
        %outside the bounds, go back to the beginning of the while loop to
        %search for a new model
        if all(x2>=xbnds(:,1)) && all(x2<=xbnds(:,2)), flag = 1; end
    end
    
    %Evaluate forward model for the proposed model x2 that is within bounds
    dhat = dhatFunc(x2);

    %Evaluate probability of the model given the data using Bayesian
    %approach in log space: 
    %log P(x2|d) = log P(x2) + log L(d|x2) 
    Px2   = PriorFunc(x2);
    Ld_x2 = LikeFunc(dhat,x2);
    Px2_d = Px2 + Ld_x2;

    %compare posterior probability with previous model
    P_accept = Px2_d-Px1_d;
    
    %some random number on the interval [0,1]
    u = log(rand(1))*anneal.temp;   
    
    %Analyze the acceptance criterion by taking the ratio of the
    %probability of the proposed model to the previous model and comparing
    %this probability to the random number u between 0 and 1.
    if P_accept >= u         
        %i.e. accept model, so store the proposed model and its likelihood
        x_keep(:,i) = x2;  
        P_keep(i) = Px2_d;
        
        %assign the accepted model for the next comparison
        x1 = x2;            
        Px1_d = Px2_d;
        count = count+1;
        
    else
        %reject this model
        x_keep(:,i) = x1;
        P_keep(i) = Px1_d;
    end
    if i>=anneal.burnin && Px2_d>=Pbest % found new best fit model
        xbest = x2;
        Pbest = Px2_d;
        ibest = i;
    end

    anneal.hist(i) = anneal.temp;

    if i==anneal.burnin; count = 0; end

    if i>anneal.burnin && i<Niter-anneal.refine
        AcceptRatio = 100*count/(i-anneal.burnin);
    else
        AcceptRatio = 100*count/(i);
    end

    if mod(i,Nprint)==0
        figure(100); clf;
        subplot(2,1,1)
        semilogy(1:i,-P_keep(1:i)); axis tight; box on; hold on; 
        ylimits = ylim;
        plot(anneal.burnin*ones(1,2), ylimits, 'Color', 0.8*[1 1 1]);
        if anneal.burnin>1; text(anneal.burnin, mean(ylimits), 'burn in', 'FontSize', 16, 'Rotation', 90, 'VerticalAlignment','bottom','HorizontalAlignment','center'); end
        plot(anneal.ilvl.*ones(2,1), ylimits.', 'Color', 0.8*[1 1 1]);
        plot((Niter-anneal.refine)*ones(1,2), ylimits, 'Color', 0.8*[1 1 1]);
        if anneal.refine>1; text((Niter-anneal.refine), mean(ylimits), 'refine', 'FontSize', 16, 'Rotation', 90, 'VerticalAlignment','top','HorizontalAlignment','center'); end
        plot(Niter*ones(1,2), ylimits, 'Color', 0.8*[1 1 1]);
        if i>=(anneal.burnin)
            plot(ibest*ones(1,2), ylimits, 'Color', 0.8*[1 0.1 0.2]); 
            text(ibest, mean(ylimits), 'best fit', 'FontSize', 16, 'Rotation', 90, 'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
        hold off
        set(gca,'ydir','reverse');
        yticklabs = get(gca,'YTickLabel');
        set(gca,'YTickLabel',strcat('-',yticklabs));
        xlabel('Iteration number','FontSize', 16); ylabel('log Likelihood', 'FontSize', 16);
        title(['Acceptance ratio = ' num2str(AcceptRatio,4)],'FontSize', 16);

        subplot(2,1,2)
        plot(1:i,anneal.hist(1:i)); axis tight; box on; hold on; 
        ylimits = [0,1];
        plot(anneal.burnin*ones(1,2), ylimits, 'Color', 0.8*[1 1 1]);
        if anneal.burnin>1; text(anneal.burnin, mean(ylimits), 'burn in', 'FontSize', 16, 'Rotation', 90, 'VerticalAlignment','bottom','HorizontalAlignment','center'); end
        plot(anneal.ilvl.*ones(2,1), ylimits.', 'Color', 0.8*[1 1 1]);
        plot((Niter-anneal.refine)*ones(1,2), ylimits, 'Color', 0.8*[1 1 1]);
        if anneal.refine>1; text((Niter-anneal.refine), mean(ylimits), 'refine', 'FontSize', 16, 'Rotation', 90, 'VerticalAlignment','top','HorizontalAlignment','center'); end
        plot(Niter*ones(1,2), ylimits, 'Color', 0.8*[1 1 1]);
        hold off;
        xlabel('Iteration number','FontSize', 16); ylabel('Annealing factor', 'FontSize', 16);
        drawnow;
    end
end

progressbar('_end')

x_keep = x_keep';

fprintf('\n\n Acceptance ratio = %.2f percent.\n\n\n', count/length(P_keep)*100);
end
