% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_default')

% test decreasing time step
RTOL = [1e-3,1e-6,1e-9];

for rtol = RTOL

    % set run parameters
    runID    =  'bnchm_cnsv';        % run identifier
    nop      =  20;                  % output frame plotted/saved every 'nop' time steps
    plot_op  =  1;                   % switch on to live plot of results
    plot_cv  =  1;                   % switch on to live plot iterative convergence

    % set model domain parameters
    D        =  10;                  % chamber depth [m]
    L        =  5;                   % chamber width [m]
    N        =  120;                 % number of grid points in z-direction (incl. 2 ghosts)
    h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

    % set model timing parameters
    Nt       =  nop;                 % number of time steps to take
    dt       =  1;                % set initial time step
    dtmax    =  1;

    % set initial thermo-chemical state
    smth     =  15;
    T0       =  1100;                % temperature top  layer [deg C]
    T1       =  T0;                  % temperature base layer [deg C]
    c0       =  [0.10  0.19  0.34  0.27  0.06  0.04  0.05];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
    c1       =  c0;                  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
    dcr      =  [1,1,1,-1,-1,-1,0]*1e-4;  % amplitude of random noise [wt]
    dcg      =  [-1,-1,-1,1,1,1,0]*1e-3;  % amplitude of centred gaussian [wt]
    dTg      =  5;
    dTr      =  0.1;

    % closed boundaries for gas flux
    fin      =  0;
    fout     =  0;

    calID    =  'MORB_hi';              % phase diagram calibration

    % set numerical model parameters
    TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
    ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
    CFL      =  100;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
    atol     =  rtol/1e6;            % outer its absolute tolerance
    maxit    =  100;                 % maximum outer its
    alpha    =  0.75;                % iterative step size parameter
    beta     =  0.00;                % iterative damping parameter
    gamma    =  0;

    % create output directory
    if ~isfolder([outdir,'/',runID])
        mkdir([outdir,'/',runID]);
    end

    % run code
    run('../src/main')

    % plot convergence
    ES = norm(diff(hist.ES(Nt/2:Nt  ))/dt,'fro')./sqrt(length(diff(hist.ES(Nt/2:Nt))         ));
    EB = norm(diff(hist.EB(Nt/2:Nt  ))/dt,'fro')./sqrt(length(diff(hist.EB(Nt/2:Nt))         ));
    % EM = norm(diff(hist.EM(Nt/2:Nt  ))/dt,'fro')./sqrt(length(diff(hist.EM(Nt/2:Nt))         ));
    % EX = norm(diff(hist.EX(Nt/2:Nt  ))/dt,'fro')./sqrt(length(diff(hist.EX(Nt/2:Nt))         ));
    % EF = norm(diff(hist.EF(Nt/2:Nt  ))/dt,'fro')./sqrt(length(diff(hist.EF(Nt/2:Nt))         ));
    EC = norm(diff(hist.EC(Nt/2:Nt,:))/dt,'fro')./sqrt(length(diff(hist.EC(Nt/2:Nt))*cal.ncmp));
    ET = norm(diff(hist.ET(Nt/2:Nt,:))/dt,'fro')./sqrt(length(diff(hist.ET(Nt/2:Nt))*cal.ntrc));

    fh20 = figure(20);
    subplot(4,1,1);
    semilogy((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,abs(diff(hist.ES(Nt/2:Nt)))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,2);
    semilogy((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,abs(diff(hist.EB(Nt/2:Nt)))/dt,'-',LW{:}); hold on; axis tight; box on; hold on
    % semilogy((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,abs(diff(hist.EM(Nt/2:Nt)))/dt,'-',LW{:}); hold on; axis tight; box on;
    % semilogy((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,abs(diff(hist.EX(Nt/2:Nt)))/dt,'-',LW{:}); hold on; axis tight; box on;
    % semilogy((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,abs(diff(hist.EF(Nt/2:Nt)))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_{\bar{\rho}}, E_{F^i}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,3);
    semilogy((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,abs(diff(hist.EC(Nt/2:Nt,:)))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_{C_j}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,4);
    semilogy((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,abs(diff(hist.ET(Nt/2:Nt,:)))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_{\Theta_k}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    xlabel('Time [s]',TX{:},FS{:});

    fh21 = figure(21);
    p1 = loglog(rtol,ES,'kd','MarkerSize',8,'LineWidth',2); hold on; box on;
    p2 = loglog(rtol,EB,'rs','MarkerSize',8,'LineWidth',2);
    % p3 = loglog(rtol,EM,'rs','MarkerSize',8,'LineWidth',2);
    % p4 = loglog(rtol,EX,'bs','MarkerSize',8,'LineWidth',2);
    % p5 = loglog(rtol,EF,'cs','MarkerSize',8,'LineWidth',2);
    p6 = loglog(rtol,EC,'go','MarkerSize',8,'LineWidth',2);
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('rel. nonlinear tolerance [1]','Interpreter','latex','FontSize',16)
    ylabel('rel. numerical error rate [1/s]','Interpreter','latex','FontSize',16)
    title('Global conservation in time','Interpreter','latex','FontSize',20)

    if rtol == RTOL(end)
        p7 = loglog(RTOL,eps.*ones(size(RTOL)),'k-' ,'LineWidth',2);  % plot trend for comparison
        legend([p1,p2,p3,p4,p5,p6,p7],{'error $S$','error $\bar{\rho}$','error $C_j$','error $\Theta_k$','machine prec.'},'Interpreter','latex','box','on','location','southeast')
        % legend([p1,p2,p3,p4,p5,p6,p7],{'error $S$','error $\bar{\rho}$','error $M$','error $X$','error $F$','error $C_j$','error $\Theta_k$','machine prec.'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

end

name = [outdir,'/',runID,'/',runID,'_',TINT];
print(fh15,name,'-dpng','-r300','-vector');