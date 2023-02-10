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
TINT     =  'bd2si';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-12;               % outer its absolute tolerance
lambda   =  0.50;                % iterative step size
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

dt    = (h/2)^2./max(kT(:)./rho(:)./cP)/2;
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

        dt = (h/2)^2./max(kT(:)./rho(:)./cP)/2;

        % solve thermo-chemical equations
        thermochem;

        % update non-linear parameters and auxiliary variables
        g0 = 10.; update;

        wx(:) = 0;  wm(:) = 0;  wf(:) = 0;
        diss(:) = 0;

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
    plot(XX(N/2,:),Sinit(N/2,:)./rhoin(N/2,:),'k',XX(N/2,:),S(N/2,:)./rho(N/2,:),'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,2)
    plot(XX(N/2,:),Tinit(N/2,:)-273.15         ,'k',XX(N/2,:),T(N/2,:)-273.15       ,'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,3)
    plot(XX(N/2,:),Cinit(N/2,:)./rhoin(N/2,:)./(1-finit(N/2,:)),'k',XX(N/2,:),C(N/2,:)./rho(N/2,:)./(1-f(N/2,:)),'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,4)
    plot(XX(N/2,:),Vinit(N/2,:)./rhoin(N/2,:)./(1-xinit(N/2,:)),'k',XX(N/2,:),V(N/2,:)./rho(N/2,:)./(1-x(N/2,:)),'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,5)
    plot(XX(N/2,:),Xinit(N/2,:)./rhoin(N/2,:),'k',XX(N/2,:),X(N/2,:)./rho(N/2,:),'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,6)
    plot(XX(N/2,:),Finit(N/2,:)./rhoin(N/2,:),'k',XX(N/2,:),F(N/2,:)./rho(N/2,:),'r','LineWidth',1.5); axis tight; box on;
    drawnow;

end

name = [opdir,'/',runID,'/',runID,'_bnchm'];
print(fh15,name,'-dpng','-r300','-vector');