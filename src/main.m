% initialise model run
init;

% physical time stepping loop
while time <= tend && step <= Nt && any(m(:)>1e-9) && ~any(cal.Tliq(:)-cal.Tsol(:)<=5) && any(T(:)-273.15-cal.Tsol(:)>1)
    
    % time step info
    timing;

    % store previous solution
    store;
    
    % reset residuals and iteration count
    resnorm  = 1;
    resnorm0 = resnorm;
    iter     = 1;
    
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

    % fractionation mode for 0D-models
    if Nx==1 && Nz==1; fractionate; end

    % record model history
    history;

    % print model diagnostics
    diagnose;

    % plot model results
    if ~mod(step,nop); output; end
    
    % increment time/step
    time = time+dt;
    step = step+1;
    if frst; frst=0; end
    
end

% save final state of model
output;

diary off
