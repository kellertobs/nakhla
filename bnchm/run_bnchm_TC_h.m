% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_default')

% set run parameters
runID    =  'bnchm_TC_h';        % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  64;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
save_op  =  0;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  4;                   % chamber depth [m]
L        =  8;                   % chamber width [m]
N        =  100;                 % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

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
TINT     =  'bd2si';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno3';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
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

NN = [20,40,80];

dt     =  L/NN(3)/16;
dtmax  =  L/NN(3)/16;

Nt     =  L/NN(1)/dt;         % number of time steps to take

for Ni = NN
    
    N     =  Ni;                % number of grid points in z-direction (incl. 2 ghosts)
    h     =  L/N;               % grid spacing (equal in both dimensions, do not set) [m]

    % run advection-diffusion benchmark
    init;

    % set velocities to constant values for lateral translation with no segregation
    W(:) = 0; Wm(:) = 0; Wx(:) = 0; Wf(:) = 0; wx(:) = 0;  wm(:) = 0;  wf(:) = 0;
    U(:) = 0; Um(:) = 1; Ux(:) = 1; Uf(:) = 1;
    P(:) = 0;

    % set periodic boundary conditions for advection
    BCA = {'periodic','periodic'};

    % set diffusion parameters to zero to isolate advection
    kT0 = 0; ks(:) = 0; kc(:) = 0; kx(:) = 0;  kf(:) = 0;

    % set parameters for non-dissipative, non-reactive flow
    diss(:) = 0;
    res_rho = 0.*rho;

    rhoin = rho; rhoout = circshift(rho,Ni/NN(1),2);
    Sin   = S;   Sout   = circshift(S  ,Ni/NN(1),2);
    Cin   = C;   Cout   = circshift(C  ,Ni/NN(1),2);
    Min   = M;   Mout   = circshift(M  ,Ni/NN(1),2);
    Xin   = X;   Xout   = circshift(X  ,Ni/NN(1),2);
    Fin   = F;   Fout   = circshift(F  ,Ni/NN(1),2);

    time  =  0;

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
            ks(:) = 0;  kc(:) = 0;  kv(:) = 0;  kx(:) = 0;  kf(:) = 0;  km(:) = 0;

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
    p1 = loglog(h,EM,'kd','MarkerSize',8,'LineWidth',2); hold on; box on;
    p2 = loglog(h,ES,'rs','MarkerSize',8,'LineWidth',2);
    p3 = loglog(h,EC,'go','MarkerSize',8,'LineWidth',2);
    p4 = loglog(h,EX,'m+','MarkerSize',8,'LineWidth',2);
    p5 = loglog(h,EF,'c^','MarkerSize',8,'LineWidth',2);
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('grid step [m]','Interpreter','latex','FontSize',16)
    ylabel('rel. numerical error [1]','Interpreter','latex','FontSize',16)
    title('Numerical convergence in space','Interpreter','latex','FontSize',20)

    if Ni == NN(1)
        p6 = loglog(L./NN,geomean([EC,ES,EX,EF]).*(NN./NN(1)).^-3,'k-','LineWidth',2);  % plot trend for comparison
        p7 = loglog(L./NN,geomean([EC,ES,EX,EF]).*(NN./NN(1)).^-5,'k--','LineWidth',2);  % plot trend for comparison
    end
    if Ni == NN(end)
        legend({'error $M$','error $S$','error $C$','error $X$','error $F$','cubic','quintic'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

    clear DD GG KP;

end

name = [opdir,'/',runID,'/',runID,'_',ADVN];
print(fh15,name,'-dpng','-r300','-vector');