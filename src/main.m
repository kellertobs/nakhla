%*****  initialise model run  *********************************************
init;

%*****  physical time stepping loop  **************************************
while time <= tend && step <= Nt &&  any(mq(:)>sqrt(eps))        ...
                                 && ~any(Tliq(:)    -Tsol(:)<=5) ...
                                 &&  any(T(:)-273.15-Tsol(:)>=5)
    
    %***  time step info
    timing;

    %***  store previous solution
    store;
    
    %***  reset residuals and iteration count
    resnorm  = 1;
    resnorm0 = resnorm;
    iter     = 1;
    if frst; alpha = alpha/2; beta = beta/2; end

    %***  non-linear iteration loop
    while resnorm/resnorm0 >= rtol/(1 + frst*10) && resnorm >= atol/(1 + frst*10) && iter <= maxit*(1 + frst)
        
        %***  solve thermo-chemical equations
        thermochem;

        %***  solve fluid-mechanics equations
        fluidmech;

        %***  update non-linear parameters and auxiliary variables
        update;

        %***  update geochemical evolution
        geochem;

        %***  report convergence
        report;

        iter = iter+1;  % increment iteration count

    end % end non-linear iterations

    %*** update phase equilibrium
    phseql;

    %***  update correlation length for convective/turbulent regularisation
    corrl;

    % renormalise sum of phase densities to bulk density
    % X = x.*rho;  M = m.*rho;  F = f.*rho;  RHO = X+M+F;

    %***  fractionation mode for 0D-models
    fractionate;

    %***  record model history
    if ~mod(step,nrh); history; end

    %***  print model diagnostics
    diagnose;

    %***  plot model results
    if ~mod(step,nop); output; end

    %***  increment time/step
    time = time+dt;
    step = step+1;
    if frst; alpha = alpha*2; beta = beta*2; frst=0; end
    
end % end time stepping

%***  save final state of model
output;

diary off
