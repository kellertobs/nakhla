% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_default')

% set run parameters
runID    =  'bnchm_TC_diff';     % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  10;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
save_op  =  0;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]

% set initial thermo-chemical state
T0       =  1140;                % temperature top layer [deg C]
dTg      =  0;                   % amplitude of centred gaussian [deg C]
c0       =  0.51;                % major component top layer [wt SiO2]
dcg      = -0.01;                % amplitude of centred gaussian [wt SiO2]
v0       =  0.04;                % volatile component top layer [wt H2O]
dvg      =  0.01;                % amplitude of centred gaussian [wt H2O]

% set model trace and isotope geochemistry parameters (must match # trace elements and isotope ratios in calibration!)
te0      =  [1,1,1,1];           % trace elements top layer [wt ppm]
dteg     =  [1,1,1,1];           % trace elements centred gaussian [wt ppm]
ir0      =  [1,1];               % isotope ratios top layer [delta]
dirg     =  [1,1];               % isotope ratios centred gaussian [delta]

finit = 0; fout = 0; Twall = nan;

% set numerical model parameters
CFL      =  1.00;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-12;               % outer its absolute tolerance
maxit    =  50;                  % maximum outer its

% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

cd ../src

N     =  100 + 2;           % number of grid points in z-direction (incl. 2 ghosts)
h     =  D/(N-2);           % grid spacing (equal in both dimensions, do not set) [m]

dt    =  1e+2;
dtmax =  1e+5;

Nt    =  100;              % number of time steps to take

% run advection-diffusion benchmark
g0 = 0.;  init;

% set velocities to zero to isolate diffusion
W(:) = 0; Wm(:) = 0; Wx(:) = 0; Wf(:) = 0; wx(:) = 0;  wm(:) = 0;  wf(:) = 0;  Vel(:) = 0;
U(:) = 0; Um(:) = 0; Ux(:) = 0; Uf(:) = 0;
P(:) = 0;

% set parameters for non-dissipative, non-reactive flow
diss(:) = 0;
tau_r   = 0;
resnorm_VP = 0;

rhoin = rho;
Tinit = T;
xinit = x;
finit = f;
Sinit = S;
Cinit = C;
Vinit = V;
Xinit = X;
Finit = F;

dt    =  1e+2;
dtmax =  1e+5;
time  =  0;

% physical time stepping loop
while time <= tend && step <= Nt

    fprintf(1,'\n\n*****  step %d;  dt = %4.4e;  time = %4.4e [hr]\n\n',step,dt./3600,time./3600);
    TTtime  = tic;
    EQtime  = 0;
    FMtime  = 0;
    TCtime  = 0;
    UDtime  = 0;

    if step==1
        alpha1 = 1; alpha2 = 1; alpha3 = 0;
        beta1  = 1; beta2  = 0; beta3  = 0;
    elseif step==2
        alpha1 = 1;   alpha2 = 1;   alpha3 = 0;
        beta1  = 1/2; beta2  = 1/2; beta3  = 0;
    else
        alpha1 = 3/2; alpha2 = 4/2; alpha3 = -1/2;
        beta1  = 3/4; beta2  = 2/4; beta3  = -1/4;
    end

    % store previous solution
    Soo = So; So = S;
    Coo = Co; Co = C;
    Voo = Vo; Vo = V;
    Xoo = Xo; Xo = X;
    Foo = Fo; Fo = F;
    TEoo = TEo; TEo = TE;
    IRoo = IRo; IRo = IR;
    dSdtoo = dSdto; dSdto = dSdt;
    dCdtoo = dCdto; dCdto = dCdt;
    dVdtoo = dVdto; dVdto = dVdt;
    dXdtoo = dXdto; dXdto = dXdt;
    dFdtoo = dFdto; dFdto = dFdt;
    dTEdtoo = dTEdto; dTEdto = dTEdt;
    dIRdtoo = dIRdto; dIRdto = dIRdt;
    rhooo = rhoo; rhoo = rho;
    Div_rhoVoo = Div_rhoVo; Div_rhoVo = Div_rhoV;
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

        % update non-linear parameters and auxiliary variables
        g0 = 10.; update;

        wx(:) = 0;  wm(:) = 0;  wf(:) = 0;
        diss(:) = 0;
        dt = (h/2)^2./max(kT(:)./rho(:)./cP)/2;

        % solve fluid-mechanics equations
        g0 = 0.; fluidmech;

        % update geochemical evolution
        geochem;

        % report convergence
        report;

        iter = iter+1;
    end

    % update history
    history;

    % plot results
    if ~mod(step,nop); output; end

    % increment time/step
    time = time+dt;
    step = step+1;

    figure(100); clf;
    subplot(3,2,1)
    plot(XX(N/2,inx),Sinit(N/2,inx)./rhoin(N/2,inx),'k',XX(N/2,inx),S(N/2,inx)./rho(N/2,inx),'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,2)
    plot(XX(N/2,inx),Tinit(N/2,inx)-273.15         ,'k',XX(N/2,inx),T(N/2,inx)-273.15       ,'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,3)
    plot(XX(N/2,inx),Cinit(N/2,inx)./rhoin(N/2,inx)./(1-finit(N/2,inx)),'k',XX(N/2,inx),C(N/2,inx)./rho(N/2,inx)./(1-f(N/2,inx)),'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,4)
    plot(XX(N/2,inx),Vinit(N/2,inx)./rhoin(N/2,inx)./(1-xinit(N/2,inx)),'k',XX(N/2,inx),V(N/2,inx)./rho(N/2,inx)./(1-x(N/2,inx)),'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,5)
    plot(XX(N/2,inx),Xinit(N/2,inx)./rhoin(N/2,inx),'k',XX(N/2,inx),X(N/2,inx)./rho(N/2,inx),'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,6)
    plot(XX(N/2,inx),Finit(N/2,inx)./rhoin(N/2,inx),'k',XX(N/2,inx),F(N/2,inx)./rho(N/2,inx),'r','LineWidth',1.5); axis tight; box on;
    drawnow;

end

name = [opdir,'/',runID,'/',runID,'_bnchm'];
print(fh15,name,'-dpng','-r300','-vector');