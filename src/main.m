% initialise model run
init;
    
% physical time stepping loop
while time <= tend && step <= Nt && any(m(:)>TINY^0.5)
    
    fprintf(1,'*****  step %d;  dt = %4.4e;  time = %4.4e [hr]\n\n',step,dt./3600,time./3600);
    TTtime  = tic;
    EQtime  = 0;
    FMtime  = 0;
    TCtime  = 0;
    UDtime  = 0;

    if     strcmp(TINT,'be1im') || step==1 % first step / 1st-order backward-Euler implicit scheme
        a1 = 1; a2 = 1; a3 = 0;
        b1 = 1; b2 = 0; b3 = 0;
    elseif strcmp(TINT,'bd2im') || step==2 % second step / 2nd-order 3-point backward-difference implicit scheme
        a1 = 3/2; a2 = 4/2; a3 = -1/2;
        b1 = 1;   b2 =  0;  b3 = 0;
    elseif strcmp(TINT,'cn2si')            % other steps / 2nd-order Crank-Nicolson semi-implicit scheme
        a1 = 1;   a2 = 1;   a3 = 0;
        b1 = 1/2; b2 = 1/2; b3 = 0;
    elseif strcmp(TINT,'bd2si')            % other steps / 2nd-order 3-point backward-difference semi-implicit scheme
        a1 = 3/2; a2 = 4/2; a3 = -1/2;
        b1 = 3/4; b2 = 2/4; b3 = -1/4;
    end

    % store previous solution
    Soo = So; So = S;
    Coo = Co; Co = C;
    Voo = Vo; Vo = V;
    Xoo = Xo; Xo = X;
    Foo = Fo; Fo = F;
    Moo = Mo; Mo = M;
    rhooo = rhoo; rhoo = rho;
    TEoo = TEo; TEo = TE;
    IRoo = IRo; IRo = IR;
    dSdtoo = dSdto; dSdto = dSdt;
    dCdtoo = dCdto; dCdto = dCdt;
    dVdtoo = dVdto; dVdto = dVdt;
    dXdtoo = dXdto; dXdto = dXdt;
    dFdtoo = dFdto; dFdto = dFdt;
    dMdtoo = dMdto; dMdto = dMdt;
    drhodtoo = drhodto; drhodto = drhodt;
    dTEdtoo = dTEdto; dTEdto = dTEdt;
    dIRdtoo = dIRdto; dIRdto = dIRdt;
    Div_Vo  = Div_V;
    dto     = dt;
    
    % reset residuals and iteration count
    resnorm  = 1;
    resnorm0 = resnorm;
    iter     = 1;
    
    % non-linear iteration loop
    while resnorm/resnorm0 >= rtol && resnorm >= atol && iter <= maxit
        
        % solve thermo-chemical equations
        thermochem;

        % solve fluid-mechanics equations
        fluidmech;

        % update non-linear parameters and auxiliary variables
        update;

        % update geochemical evolution
        geochem;

        % report convergence
        report;

        iter = iter+1;
    end
    
    % record model history
    history;

    % print diagnostics
    fprintf(1,'\n         total time to solution = %3.3f sec\n\n',toc(TTtime));
    fprintf(1,'         thermo-chemical solve  = %1.3e sec\n'  ,TCtime/(iter-1));
    fprintf(1,'         phase equilibr. solve  = %1.3e sec\n'  ,EQtime/(iter-1));
    fprintf(1,'         coefficients update    = %1.3e sec\n'  ,UDtime/(iter-1));
    fprintf(1,'         fluid-mechanics solve  = %1.3e sec\n\n',FMtime/(iter-1));
    
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

% save final state of model
output;

diary off
