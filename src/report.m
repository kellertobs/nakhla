if it<=0 || resnorm>resnorm0; resnorm0 = resnorm + 1e-32; end  % reset reference residual

% report iterations
if     it >=  0  && it <  10
    fprintf(1,'    ---  it =    %d;   abs res = %4.4e;   rel res = %4.4e \n',it,resnorm,resnorm/resnorm0);
elseif it >= 10  && it < 100
    fprintf(1,'    ---  it =   %d;   abs res = %4.4e;   rel res = %4.4e \n',it,resnorm,resnorm/resnorm0);
elseif it >= 100 && it < 1000
    fprintf(1,'    ---  it =  %d;   abs res = %4.4e;   rel res = %4.4e \n',it,resnorm,resnorm/resnorm0);
end 

% plot convergence of outer iterations
if plot_cv
    figure(100); if it==0; clf; else; hold on; end
    plot(it,log10(resnorm),'r.','MarkerSize',15); box on; axis tight;
    drawnow;
end
