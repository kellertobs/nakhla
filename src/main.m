% initialise model run
init;

% physical time stepping loop
while time <= tend && step <= M
        
    fprintf(1,'\n\n\n*****  step %d;  dt = %4.4e;  time = %4.4e [hr]\n\n',step,dt./3600,time./3600);
    tic;
    
    % store previous solution
    Ho      = H;
    Co      = C;
    Vo      = V;
    To      = T;
    co      = c;
    vo      = v;
    xo      = x;
    fo      = f;
    mo      = m;
    SImo    = SIm;
    SIxo    = SIx;
    ITo     = IT;
    CTo     = CT;
    rhoo    = rho;
    etao    = eta;
    dHdto   = dHdt;
    dCdto   = dCdt;
    dVdto   = dVdt;
    dxdto   = dxdt;
    dfdto   = dfdt;
    dITdto  = dITdt;
    dCTdto  = dCTdt;
    dSImdto = dSImdt;
    dSIxdto = dSIxdt;
    Div_rhoVo = Div_rhoV;
    Pto     = Pt;
    dto     = dt;
    
    % reset residuals and iteration count
    resnorm  = 1e3;
    resnorm0 = resnorm;
    iter     = 0;
    
    % non-linear iteration loop
    while resnorm/resnorm0 >= rtol && resnorm >= atol && iter <= maxit || iter <= 2
            
        % solve thermo-chemical equations
        thermochem;
        
        % update non-linear parameters and auxiliary variables
        update;
        
        % solve fluid-mechanics equations
        fluidmech;
        
        % update non-linear parameters and auxiliary variables
        update;
        
        % report convergence
        report;

        iter = iter+1;
    end
    
    % print diagnostics
    fprintf(1,'\n         time to solution = %4.4f sec\n\n',toc);
    
    fprintf(1,'         min T   =  %4.1f;    mean T   = %4.1f;    max T   = %4.1f;   [degC]\n' ,min(T(:)  ),mean(T(:)  ),max(T(:)  ));
    fprintf(1,'         min c   =  %1.4f;    mean c   = %1.4f;    max c   = %1.4f;   [wt]\n'   ,min(c(:)  ),mean(c(:)  ),max(c(:)  ));
    fprintf(1,'         min v   =  %1.4f;    mean v   = %1.4f;    max v   = %1.4f;   [wt]\n\n' ,min(v(:)  ),mean(v(:)  ),max(v(:)  ));
    
    fprintf(1,'         min x   =  %1.4f;    mean x   = %1.4f;    max x   = %1.4f;   [wt]\n'   ,min(x(:)  ),mean(x(:)  ),max(x(:)  ));
    fprintf(1,'         min f   =  %1.4f;    mean f   = %1.4f;    max f   = %1.4f;   [wt]\n'   ,min(f(:)  ),mean(f(:)  ),max(f(:)  ));
    fprintf(1,'         min m   =  %1.4f;    mean m   = %1.4f;    max m   = %1.4f;   [wt]\n\n' ,min(m(:)  ),mean(m(:)  ),max(m(:)  ));

    fprintf(1,'         min rho =  %4.1f;    mean rho = %4.1f;    max rho = %4.1f;   [kg/m3]\n'  ,min(rho(:)),mean(rho(:))   ,max(rho(:)));
    fprintf(1,'         min eta =  %1.2e;  mean eta = %1.2e;  max eta = %1.2e; [Pas]\n\n',min(eta(:)),geomean(eta(:)),max(eta(:)));

    fprintf(1,'         min U   = %1.4f;    mean U   = %1.4f;    max U   = %1.4f;   [m/s]\n'  ,min(U(:)  ),mean(U(:)  ),max(U(:)  ));
    fprintf(1,'         min W   = %1.4f;    mean W   = %1.4f;    max W   = %1.4f;   [m/s]\n'  ,min(-W(:) ),mean(-W(:) ),max(-W(:) ));
    fprintf(1,'         min P   = %2.4f;    mean P   = %2.4f;    max P   = %2.4f;  [kPa]\n\n',min(P(:)./1e3),mean(P(:)./1e3),max(P(:)./1e3));

    % plot results
    if ~mod(step,nop); output; end
    
    % increment time/step
    time = time+dt;
    step = step+1;
    
end

diary off
