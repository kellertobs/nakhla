% initialise model run
init;

if bnchm
    % run manufactured solution benchmark on fluid mechanics solver if specified
    mms;
    return;
end
    
% physical time stepping loop
while time <= tend && step <= M
        
    fprintf(1,'*****  step %d;  dt = %4.4e;  time = %4.4e [hr]\n\n',step,dt./3600,time./3600);
    tic;

    if step==1; theta = 1; else; theta = 0.5; end

    % store previous solution
    So      = S;
    Co      = C;
    Vo      = V;
    To      = T;
    co      = c;
    vo      = v;
    Xo      = X;
    Fo      = F;
    xo      = x;
    fo      = f;
    mo      = m;
    ITo     = IT;
    CTo     = CT;
    SIo     = SI;
    RIPo    = RIP;
    RIDo    = RID;
    rhoo    = rho;
    Div_rhoVo =  Div_rhoV;
    etao    = eta;
    dSdto   = dSdt;
    dCdto   = dCdt;
    dVdto   = dVdt;
    dXdto   = dXdt;
    dFdto   = dFdt;
    dITdto  = dITdt;
    dCTdto  = dCTdt;
    dSIdto  = dSIdt;
    dRIPdto = dRIPdt;
    dRIDdto = dRIDdt;
    Pto     = Pt;
    dto     = dt;
    
    % reset residuals and iteration count
    resnorm  = 1;
    resnorm0 = resnorm;
    iter     = 1;
    
    % non-linear iteration loop
    while resnorm/resnorm0 >= rtol && resnorm >= atol && iter <= maxit
        
        % solve thermo-chemical equations
        thermochem;
                
        % update non-linear parameters and auxiliary variables
        update;

        % solve fluid-mechanics equations
        fluidmech;
        
        % report convergence
        report;

        % update geochemical evolution
        geochem;

        iter = iter+1;
    end
    
    % record model history
    history;
    
    % print diagnostics
    fprintf(1,'\n         time to solution = %4.4f sec\n\n',toc);
    
    fprintf(1,'         min T   =  %4.1f;    mean T   = %4.1f;    max T   = %4.1f;   [degC]\n' ,min(T(:)-273.15),mean(T(:)-273.15),max(T(:)-273.15));
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
