% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_default')

% set run parameters
runID    =  'bnchm_TC_dt';       % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  80;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
save_op  =  0;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  4;                   % chamber depth [m]
L        =  8;                   % chamber width [m]
N        =  80;                 % number of grid points in z-direction (incl. 2 ghosts)
h        =  L/N;                 % grid spacing (equal in both dimensions, do not set) [m]

% set initial thermo-chemical state
T0       =  1050;                % temperature top  layer [deg C]
T1       =  T0;                  % temperature base layer [deg C]
c0       =  [0.04  0.12  0.44  0.24  0.14  0.02  0.04];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1       =  c0;                  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr      =  [1,1,1,-1,-1,-1,0]*0e-4;  % amplitude of random noise [wt]
dcg      =  [-1,-1,-1,1,1,1,0]*1e-2;  % amplitude of centred gaussian [wt]
dg_trc   =  [-1,-1,-1,1,1,1]  *1e-2;  % trace elements centred gaussian [wt ppm]

fin = 0; fout = 0; Twall = [nan,nan,nan]; periodic = 1;

% set numerical model parameters
CFL      =  1;                   % (physical) time stepping courant number (multiplies stable step) [0,1]
TINT     =  'cn2si';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
rtol     =  1e-9;                % outer its relative tolerance
atol     =  1e-15;               % outer its absolute tolerance
tau_r    =  1e32;                % disable reaction
alpha    =  1.00;                % iterative step size parameter
beta     =  0.00;                % iterative damping parameter

% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

cd ../src

DT = [h/4,h/8,h/16];

for dti = DT
    
    dt    =  dti;
    dtmax =  dti;
    Nt    =  2*h/dti;

    % initialise fields
    init;

    % set velocities to constant values for lateral translation with no segregation
    W(:) = 0; Wm(:) = 0; Wx(:) = 0; Wf(:) = 0; wx(:) = 0;  wm(:) = 0;  wf(:) = 0;
    U(:) = 0; Um(:) = 1; Ux(:) = 1; Uf(:) = 1;
    P(:) = 0;

    % set periodic boundary conditions for advection
    BCA = {'closed','periodic'};

    % set diffusion parameters to zero to isolate advection
    kT0 = 0; ks(:) = 0; kc(:) = 0; kx(:) = 0;  kf(:) = 0;

    % set parameters for non-dissipative, non-reactive flow
    diss(:) = 0;
    res_rho = 0.*rho;

    rhoin = rho; rhoout = circshift(rho,2,2);
    Sin   = S;   Sout   = circshift(S  ,2,2);
    Cin   = C;   Cout   = circshift(C  ,2,2);
    Min   = M;   Mout   = circshift(M  ,2,2);
    Xin   = X;   Xout   = circshift(X  ,2,2);
    Fin   = F;   Fout   = circshift(F  ,2,2);

    dt    = dti;
    dtmax = dti;
    time  = 0;

    output;

    % physical time stepping loop
    while time <= tend && step <= Nt

        fprintf(1,'\n\n*****  step %d;  dt = %4.4e;  time = %4.4e [hr]\n\n',step,dt./3600,time./3600);
        TTtime  = tic;
        EQtime  = 0;
        FMtime  = 0;
        TCtime  = 0;
        UDtime  = 0;

        if     strcmp(TINT,'be1im') || step==1         % first step / 1st-order backward-Euler implicit scheme
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

            % solve thermo-chemical equations
            thermochem;

            % update non-linear parameters and auxiliary variables
            update;

            wx(:) = 0;  wm(:) = 0;  wf(:) = 0;
            ks(:) = 0;  kc(:) = 0;  kx(:) = 0;  kf(:) = 0;  km(:) = 0;  kT0 = 0;
            
            [grdTx ,grdTz ] = gradient(T(icz,icx),h);
            diss = kT0./T.*(grdTz(2:end-1,2:end-1).^2 + grdTx(2:end-1,2:end-1).^2);

            % update geochemical evolution
            geochem;

            % report convergence
            report;

            iter = iter+1;
        end

        % plot results
        if ~mod(step,nop); output; end

        % increment time/step
        time = time+dt;
        step = step+1;

        figure(100); clf;
        subplot(2,1,1)
        plot(XX(ceil(N/4),:),Xout(ceil(N/4),:)./rhoout(ceil(N/4),:),'k',XX(ceil(N/4),:),X(ceil(N/4),:)./rho(ceil(N/4),:),'r','LineWidth',1.5); axis tight; box on;
        subplot(2,1,2)
        plot(XX(ceil(N/4),:),Sout(ceil(N/4),:)./rhoout(ceil(N/4),:),'k',XX(ceil(N/4),:),S(ceil(N/4),:)./rho(ceil(N/4),:),'r','LineWidth',1.5); axis tight; box on;
        drawnow;

    end


    % plot convergence
    EM = norm(rho-rhoout,'fro')./norm(rhoout,'fro');
    ES = norm(S-Sout,'fro')./norm(Sout,'fro');
    EC = norm(C-Cout,'fro')./norm(Cout,'fro');
    EX = norm(X-Xout,'fro')./norm(Xout,'fro');
    EF = norm(F-Fout,'fro')./norm(Fout,'fro');

    fh15 = figure(15);
    p1 = loglog(dt,EM,'kd','MarkerSize',8,'LineWidth',2); hold on; box on;
    p2 = loglog(dt,ES,'rs','MarkerSize',8,'LineWidth',2);
    p3 = loglog(dt,EC,'go','MarkerSize',8,'LineWidth',2);
    p4 = loglog(dt,EX,'m+','MarkerSize',8,'LineWidth',2);
    p5 = loglog(dt,EF,'c^','MarkerSize',8,'LineWidth',2);
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('time step [s]','Interpreter','latex','FontSize',16)
    ylabel('rel. numerical error [1]','Interpreter','latex','FontSize',16)
    title('Numerical convergence in time','Interpreter','latex','FontSize',20)

    if dt == DT(1)
        p6 = loglog(DT,geomean([EM,ES,EC,EX,EF]).*(DT./DT(1)).^1,'k--','LineWidth',2);  % plot trend for comparison
        p7 = loglog(DT,geomean([EM,ES,EC,EX,EF]).*(DT./DT(1)).^2,'k-' ,'LineWidth',2);
    end
    if dt == DT(end)
        legend({'error $M$','error $S$','error $C$','error $X$','error $F$','quadratic'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

end

name = [opdir,'/',runID,'/',runID,'_',TINT];
print(fh15,name,'-dpng','-r300','-vector');