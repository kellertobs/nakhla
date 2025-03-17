% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_default')

% set run parameters
runID    =  'bnchm_TC_dt';       % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  10;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
plot_cv  =  0;                   % switch on to live plot iterative convergence
save_op  =  0;

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  50;                  % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  nop;                 % number of time steps to take
dt       =  1;                   % set initial time step

% set initial thermo-chemical state
smth     =  15;
T0       =  1150;                % temperature top  layer [deg C]
T1       =  T0;                  % temperature base layer [deg C]
c0       =  [11  17  35  31  3  3  5]/100;  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1       =  c0;                  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr      =  [1,1,1,-1,-1,-1,0]*0e-3;  % amplitude of random noise [wt]
dcg      =  [-1,-1,-1,1,1,1,0]*1e-2;  % amplitude of centred gaussian [wt]
dTg      =  5;
dTr      =  0.0;
dr_trc   =  [1,1,1,-1,-1,-1].*0e-3;
dg_trc   =  [-1,-1,-1,1,1,1].*1e-2;

% closed boundaries for gas flux
periodic =  1;
bndmode  =  0;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
Ptop     =  2.0e8;               % top pressure [Pa]
fin      =  0;
fout     =  0;
tau_r    =  1e32;
calID    =  'DEMO';              % phase diagram calibration

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  1;                   % (physical) time stepping courant number (multiplies stable step) [0,1]
atol     =  1e-12;               % outer its absolute tolerance
rtol     =  atol/1e6;            % outer its absolute tolerance
maxit    =  100;                 % maximum outer its
alpha    =  0.75;                % iterative step size parameter

% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

cd ../src

DT    = [h/2,h/4,h/8];  % time step sizes to test relative to grid step
nshft = 4;              % number of grid steps target is shifted from initial

for dti = DT
    
    dt    =  dti;
    dtmax =  dti;
    Nt    =  nshft*h/dti;

    % initialise fields
    init;

    % set velocities to constant values for lateral translation with no segregation
    W(:) = 0; Wm(:) = 0; Wx(:) = 0; Wf(:) = 0; wx(:) = 0;  wm(:) = 0;  wf(:) = 0;
    U(:) = 0; Um(:) = 1; Ux(:) = 1; Uf(:) = 1;
    P(:) = 0;

    % set diffusion parameters to zero to isolate advection
    kT(:) = 0;  ks(:) = 0;  kc(:) = 0;  km(:) = 0;  kx(:) = 0;  kf(:) = 0;  diss(:) = 0;

    % set parameters for non-dissipative, non-reactive flow
    diss(:) = 0;
    res_rho = 0.*rho;

    Sin   = S;   Sout   = circshift(S  ,nshft,2);
    rhoin = rho; rhoout = circshift(rho,nshft,2);
    Min   = M;   Mout   = circshift(M  ,nshft,2);
    Xin   = X;   Xout   = circshift(X  ,nshft,2);
    Fin   = F;   Fout   = circshift(F  ,nshft,2);
    Cin   = C;   Cout   = circshift(C  ,nshft,2);
    TRCin = TRC; TRCout = circshift(TRC,nshft,2);

    dt    = dti;
    dtmax = dti;
    time  = 0;

    output;

    % physical time stepping loop
    while time <= tend && step <= Nt

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
        while resnorm/resnorm0 >= rtol && resnorm >= atol && iter <= maxit

            % solve thermo-chemical equations
            thermochem;

            % update non-linear parameters and auxiliary variables
            update;

            kT(:) = 0;  ks(:) = 0;  kc(:) = 0;  km(:) = 0;  kx(:) = 0;  kf(:) = 0;  diss(:) = 0;

            % update geochemical evolution
            geochem;

            % report convergence
            report;

            iter = iter+1;
        end

        % X = X./RHO.*rho;  M = M./RHO.*rho;  F = F./RHO.*rho;  RHO = X+M+F;

        % print model diagnostics
        diagnose;

        % plot model results
        if ~mod(step,nop); output; end

        % increment time/step
        time = time+dt;
        step = step+1;
        if frst; alpha = alpha*2; beta = beta*2; frst=0; end

        figure(100); clf;
        subplot(2,1,1)
        plot(XX(ceil(Nz/2),:),Xout(ceil(Nz/2),:)./rhoout(ceil(Nz/2),:),'k',XX(ceil(Nz/4),:),X(ceil(Nz/2),:)./rho(ceil(Nz/2),:),'r','LineWidth',1.5); axis tight; box on;
        subplot(2,1,2)
        plot(XX(ceil(Nz/2),:),Sout(ceil(Nz/2),:)./rhoout(ceil(Nz/2),:),'k',XX(ceil(Nz/2),:),S(ceil(Nz/2),:)./rho(ceil(Nz/2),:),'r','LineWidth',1.5); axis tight; box on;
        drawnow;

    end


    % plot convergence
    ES = norm(S-Sout,'fro')./norm(Sout,'fro');
    EB = norm(rho-rhoout,'fro')./norm(rhoout,'fro');
    EM = norm(M-Mout,'fro')./norm(Mout,'fro');
    EX = norm(X-Xout,'fro')./norm(Xout,'fro');
    EF = norm(F-Fout,'fro')./norm(Fout,'fro');
    EC = norm(C-Cout,'fro')./norm(Cout,'fro');
    ET = norm(TRC-TRCout,'fro')./norm(TRCout,'fro');

    clist = [colororder;[0 0 0]];

    fh15 = figure(15);
    p1 = loglog(dt,ES,'+','Color',clist(1,:),'MarkerSize',10,'LineWidth',2); hold on; box on;
    % p2 = loglog(dt,EB,'s','Color',clist(2,:),'MarkerSize',10,'LineWidth',2);
    p3 = loglog(dt,EM,'o','Color',clist(3,:),'MarkerSize',10,'LineWidth',2);
    p4 = loglog(dt,EX,'d','Color',clist(4,:),'MarkerSize',10,'LineWidth',2);
    p5 = loglog(dt,EF,'*','Color',clist(5,:),'MarkerSize',10,'LineWidth',2);
    p6 = loglog(dt,EC,'^','Color',clist(6,:),'MarkerSize',10,'LineWidth',2);
    p7 = loglog(dt,ET,'v','Color',clist(7,:),'MarkerSize',10,'LineWidth',2);
    xlabel('time step [s]','Interpreter','latex','FontSize',16)
    ylabel('rel. numerical error [1]','Interpreter','latex','FontSize',16)
    title('Numerical convergence in time','Interpreter','latex','FontSize',20)

    if dt == DT(1)
        p8 = loglog(DT,geomean([ES,EM,EX,EF,EC,ET]).*(DT./DT(1)).^1,'k--','LineWidth',2);  % plot trend for comparison
        p9 = loglog(DT,geomean([ES,EM,EX,EF,EC,ET]).*(DT./DT(1)).^2,'k-' ,'LineWidth',2);
    end
    if dt == DT(end)
        legend({'error $S$','error $M$','error $X$','error $F$','error $C_j$','error $\Theta_k$','linear','quadratic'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

end

name = [opdir,'/',runID,'/',runID,'_',TINT];
print(fh15,name,'-dpng','-r300','-vector');