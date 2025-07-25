% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_default')

% set run parameters
runID    =  'bnchm_TC_h';        % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
plot_cv  =  0;                   % switch on to live plot iterative convergence
save_op  =  0;

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                   % chamber width [m]
N        =  200;                 % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

% set initial thermo-chemical state
smth     =  15;
T0       =  1100;                % temperature top  layer [deg C]
T1       =  T0;                  % temperature base layer [deg C]
c0       =  [10  19  32  30  5  4  6]/100;  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
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
alpha    =  0.50;                % iterative step size parameter

% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

cd ../src

NN    = 25*[1,2,4];
nshft = 1;

dt     =  D/NN(3)/100;
dtmax  =  D/NN(3)/100;

Nt     =  nshft*D/NN(1)/dt;           % number of time steps to take

for Ni = NN
    
    N     =  Ni;                % number of grid points in z-direction
    h     =  D/N;               % grid spacing

    % initialise fields
    init;

    % set velocities to constant values for lateral translation with no segregation
    W(:) = 0; Wm(:) = 0; Wx(:) = 0; Wf(:) = 0; wx(:) = 0;  wm(:) = 0;  wf(:) = 0;
    U(:) = 0; Um(:) = 1; Ux(:) = 1; Uf(:) = 1;
    P(:) = 0;

    % set diffusion parameters to zero to isolate advection
    kT(:) = 0;  ks(:) = 0;  kc(:) = 0;  km(:) = 0;  kx(:) = 0;  kf(:) = 0;  diss(:) = 0;  Delta_cnv = 0;

    % set parameters for non-dissipative, non-reactive flow
    diss(:) = 0;
    res_rho = 0.*rho;

    Sin   = S;   Sout   = circshift(S  ,Ni/NN(1)*nshft,2);
    rhoin = rho; rhoout = circshift(rho,Ni/NN(1)*nshft,2);
    Min   = M;   Mout   = circshift(M  ,Ni/NN(1)*nshft,2);
    Xin   = X;   Xout   = circshift(X  ,Ni/NN(1)*nshft,2);
    Fin   = F;   Fout   = circshift(F  ,Ni/NN(1)*nshft,2);
    Cin   = C;   Cout   = circshift(C  ,Ni/NN(1)*nshft,2);
    TRCin = TRC; TRCout = circshift(TRC,Ni/NN(1)*nshft,2);

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

            kT(:) = 0;  ks(:) = 0;  kc(:) = 0;  km(:) = 0;  kx(:) = 0;  kf(:) = 0;  diss(:) = 0;  Delta_cnv = 0;

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

        if ~mod(step,100)
        figure(100); clf;
        subplot(2,1,1)
        plot(XX(ceil(Nz/2),:),Xout(ceil(Nz/2),:)./rhoout(ceil(Nz/2),:),'k',XX(ceil(Nz/4),:),X(ceil(Nz/2),:)./rho(ceil(Nz/2),:),'r','LineWidth',1.5); axis tight; box on;
        subplot(2,1,2)
        plot(XX(ceil(Nz/2),:),Sout(ceil(Nz/2),:)./rhoout(ceil(Nz/2),:),'k',XX(ceil(Nz/2),:),S(ceil(Nz/2),:)./rho(ceil(Nz/2),:),'r','LineWidth',1.5); axis tight; box on;
        drawnow;
        end
    end


    % plot convergence
    ES = norm(S-Sout,'fro')./norm(Sout,'fro');
    EB = norm(RHO-rhoout,'fro')./norm(rhoout,'fro');
    EM = norm(M-Mout,'fro')./norm(Mout,'fro');
    EX = norm(X-Xout,'fro')./norm(Xout,'fro');
    EF = norm(F-Fout,'fro')./norm(Fout,'fro');
    EC = norm(C-Cout,'fro')./norm(Cout,'fro');
    ET = norm(TRC-TRCout,'fro')./norm(TRCout,'fro');

    clist = [colororder;[0 0 0]];

    fh15 = figure(15);
    p1 = loglog(h,ES,'+','Color',clist(1,:),'MarkerSize',10,'LineWidth',2); hold on; box on;
    p2 = loglog(h,EB,'s','Color',clist(2,:),'MarkerSize',10,'LineWidth',2);
    p3 = loglog(h,EM,'o','Color',clist(3,:),'MarkerSize',10,'LineWidth',2);
    p4 = loglog(h,EX,'d','Color',clist(4,:),'MarkerSize',10,'LineWidth',2);
    p5 = loglog(h,EF,'*','Color',clist(5,:),'MarkerSize',10,'LineWidth',2);
    p6 = loglog(h,EC,'^','Color',clist(6,:),'MarkerSize',10,'LineWidth',2);
    p7 = loglog(h,ET,'v','Color',clist(7,:),'MarkerSize',10,'LineWidth',2);
    xlabel('grid step [m]','Interpreter','latex','FontSize',16)
    ylabel('rel. numerical error [1]','Interpreter','latex','FontSize',16)
    title('Numerical convergence in space','Interpreter','latex','FontSize',20)

    if Ni == NN(1)
        p8  = loglog(D./NN,geomean([ES,EB,EM,EX,EF,EC,ET]).*(NN./NN(1)).^-1,'k-.','LineWidth',2);  % plot trend for comparison
        p9  = loglog(D./NN,geomean([ES,EB,EM,EX,EF,EC,ET]).*(NN./NN(1)).^-3,'k--','LineWidth',2);  % plot trend for comparison
        p10 = loglog(D./NN,geomean([ES,EB,EM,EX,EF,EC,ET]).*(NN./NN(1)).^-5,'k-','LineWidth',2);  % plot trend for comparison
    end
    if Ni == NN(end)
        legend({'error $S$','error $\bar{\rho}$','error $M$','error $X$','error $F$','error $C_j$','error $\Theta_k$','linear','cubic','quintic'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

    clear DD GG KP;

end

name = [opdir,'/',runID,'/',runID,'_',ADVN];
print(fh15,name,'-dpng','-r300','-vector');