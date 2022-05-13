% get combined residual norm
resnorm = resnorm_VP + resnorm_TC;
if iter==1 || resnorm>resnorm0; resnorm0 = resnorm + 1e-32; end  % reset reference residual

% check for solver divergence or failing
if isnan(resnorm); error('!!! Solver failed with NaN: end run !!!'); end

% report iterations
if     iter >=  0  && iter <  10
    fprintf(1,'    ---  iter =    %d;   abs res = %4.4e;   rel res = %4.4e \n',iter,resnorm,resnorm/resnorm0);
elseif iter >= 10  && iter < 100
    fprintf(1,'    ---  iter =   %d;   abs res = %4.4e;   rel res = %4.4e \n',iter,resnorm,resnorm/resnorm0);
elseif iter >= 100 && iter < 1000
    fprintf(1,'    ---  iter =  %d;   abs res = %4.4e;   rel res = %4.4e \n',iter,resnorm,resnorm/resnorm0);
end 

% plot convergence of outer iterations
if plot_cv
    figure(100); if iter==1; clf; else; hold on; end
    plot(iter,log10(resnorm_TC),'b.',iter,log10(resnorm_VP),'r.',iter,log10(resnorm),'k.','MarkerSize',15,'LineWidth',1.5); box on; axis tight;
    drawnow;
end
