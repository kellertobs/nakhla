% initialise model run
init;

% physical time stepping loop
while time <= tend && step <= Nt &&  any(mq(:)>eps^0.5) ...
                                 && ~any(cal.Tliq(:)-cal.Tsol(:)<=3) && any(T(:)-273.15-cal.Tsol(:)>=3)
    
    % time step info
    timing;

    % store previous solution
    store;
    
    % reset residuals and iteration count
    resnorm  = 1;
    resnorm0 = resnorm;
    iter     = 1;
    if frst; alpha = alpha/2; beta = beta/2; end

    % non-linear iteration loop
    while resnorm/resnorm0 >= rtol/(1 + frst*10) && resnorm >= atol/(1 + frst*10) && iter <= maxit*(1 + frst)
        
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

    % renormalise sum of phase, component densities to bulk density
    X = X./RHO.*rho;  M = M./RHO.*rho;  F = F./RHO.*rho;  RHO = X+M+F;
    C = C./sum(C,3).*rho;

    % fractionation mode for 0D-models
    if Nx==1 && Nz==1
        Ptop = Ptop + (T-To).*dPdT;
        fractionate; 
    end

    % record model history
    history;

    % print model diagnostics
    diagnose;

    % plot model results
    if ~mod(step,nop); output; end
    
    % increment time/step
    time = time+dt;
    step = step+1;
    if frst; alpha = alpha*2; beta = beta*2; frst=0; end
    
end

% save final state of model
output;

diary off
