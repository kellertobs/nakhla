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

% set initial thermo-chemical state
T0       =  1050;                % temperature top  layer [deg C]
T1       =  T0;                  % temperature base layer [deg C]
dTg      =  0;                   % amplitude of centred gaussian [deg C]
c0       =  [0.04,0.12,0.44,0.24,0.14,0.02,0.04];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1       =  c0;                  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr      =  [0,0,0,0,0,0,0];     % amplitude of random noise [wt]
dcg      =  [-1,-1,-1,1,1,1,0]*1e-2;  % amplitude of centred gaussian [wt]

fin = 0; fout = 0; Twall = [nan,nan,nan];

% set numerical model parameters
CFL      =  1.00;                % (physical) time stepping courant number (multiplies stable step) [0,1]
TINT     =  'bd2si';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-12;               % outer its absolute tolerance
maxit    =  50;                  % maximum outer its
alpha    =  0.75;
beta     =  0.10;

% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

cd ../src

N     =  100 ;           % number of grid points in z-direction (incl. 2 ghosts)
h     =  D/N;            % grid spacing (equal in both dimensions, do not set) [m]

Nt    =  100;            % number of time steps to take

% run advection-diffusion benchmark
init;

rhoin = rho;
Tinit = T;
xinit = x;
finit = f;
Sinit = S;
Cinit = C;
Xinit = X;
Minit = M;
Finit = F;

% dt    = (h/2)^2./max([kc(:)./rho(:);(kT0+ks(:).*T(:))./rho(:)./cP])/16;
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

    if     strcmp(TINT,'be1im') || step==1 || frst % first step / 1st-order backward-Euler implicit scheme
        a1 = 1; a2 = 1; a3 = 0;
        b1 = 1; b2 = 0; b3 = 0;
    elseif strcmp(TINT,'bd2im') || step==2         % second step / 2nd-order 3-point backward-difference implicit scheme
        a1 = 3/2; a2 = 4/2; a3 = -1/2;
        b1 = 1;   b2 =  0;  b3 = 0;
    elseif strcmp(TINT,'cn2si')                    % other steps / 2nd-order Crank-Nicolson semi-implicit scheme
        a1 = 1;   a2 = 1;   a3 = 0;
        b1 = 1/2; b2 = 1/2; b3 = 0;
    elseif strcmp(TINT,'bd2si')                    % other steps / 2nd-order 3-point backward-difference semi-implicit scheme
        a1 = 3/2; a2 = 4/2; a3 = -1/2;
        b1 = 3/4; b2 = 2/4; b3 = -1/4;
    end

    % store previous solution
    Soo = So; So = S;
    Coo = Co; Co = C;
    Xoo = Xo; Xo = X;
    Foo = Fo; Fo = F;
    Moo = Mo; Mo = M;
    rhooo = rhoo; rhoo = rho;
    TRCoo = TRCo; TRCo = TRC;
    dSdtoo = dSdto; dSdto = dSdt;
    dCdtoo = dCdto; dCdto = dCdt;
    dXdtoo = dXdto; dXdto = dXdt;
    dFdtoo = dFdto; dFdto = dFdt;
    dMdtoo = dMdto; dMdto = dMdt;
    drhodtoo = drhodto; drhodto = drhodt;
    dTRCdtoo = dTRCdto; dTRCdto = dTRCdt;
    Div_Vo  = Div_V;
    rhoWoo  = rhoWo; rhoWo = rhofz.*W(:,2:end-1);
    rhoUoo  = rhoUo; rhoUo = rhofx.*U(2:end-1,:);
    Pchmboo = Pchmbo; Pchmbo = Pchmb;
    dPchmbdtoo = dPchmbdto; dPchmbdto = dPchmbdt;
    dto     = dt;

    % reset residuals and iteration count
    resnorm  = 1;
    resnorm0 = resnorm;
    iter     = 1;

    % non-linear iteration loop
    while resnorm/resnorm0 >= rtol && resnorm >= atol && iter <= maxit

        % dt = (h/2)^2./max([kc(:)./rho(:);(kT0+ks(:).*T(:))./rho(:)./cP])/16;

        W(:) = 0; Wm(:) = 0; Wx(:) = 0; Wf(:) = 0; wx(:) = 0;  wm(:) = 0;  wf(:) = 0;  Vel(:) = 0;
        U(:) = 0; Um(:) = 0; Ux(:) = 0; Uf(:) = 0;
        P(:) = 0;
        diss(:) = 0;

        % solve thermo-chemical equations
        thermochem;

        % solve fluid-mechanics equations
        fluidmech; update;

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
    plot(XX(N/2,:),Cinit(N/2,:,3)./rhoin(N/2,:)./(1-finit(N/2,:)),'k',XX(N/2,:),C(N/2,:,3)./rho(N/2,:)./(1-f(N/2,:)),'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,4)
    plot(XX(N/2,:),Cinit(N/2,:,4)./rhoin(N/2,:)./(1-xinit(N/2,:)),'k',XX(N/2,:),C(N/2,:,4)./rho(N/2,:)./(1-x(N/2,:)),'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,5)
    plot(XX(N/2,:),Xinit(N/2,:)./rhoin(N/2,:),'k',XX(N/2,:),X(N/2,:)./rho(N/2,:),'r','LineWidth',1.5); axis tight; box on;
    subplot(3,2,6)
    plot(XX(N/2,:),Finit(N/2,:)./rhoin(N/2,:),'k',XX(N/2,:),F(N/2,:)./rho(N/2,:),'r','LineWidth',1.5); axis tight; box on;
    drawnow;

end

name = [opdir,'/',runID,'/',runID,'_bnchm'];
print(fh15,name,'-dpng','-r300','-vector');