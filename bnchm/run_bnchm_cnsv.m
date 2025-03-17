% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_default')

% test decreasing time step
ATOL = [1e-6,1e-9,1e-12];

for atol = ATOL

    % set run parameters
    runID    =  'bnchm_cnsv';        % run identifier
    nop      =  20;                  % output frame plotted/saved every 'nop' time steps
    plot_op  =  1;                   % switch on to live plot of results
    plot_cv  =  0;                   % switch on to live plot iterative convergence
    save_op  =  0;

    % set model domain parameters
    D        =  10;                  % chamber depth [m]
    L        =  10;                  % chamber width [m]
    N        =  100;                 % number of grid points in z-direction (incl. 2 ghosts)
    h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

    % set model timing parameters
    Nt       =  nop;                 % number of time steps to take
    dt       =  1;                   % set initial time step

    % set initial thermo-chemical state
    smth     =  15;
    T0       =  1200;                % temperature top  layer [deg C]
    T1       =  T0;                  % temperature base layer [deg C]
    c0       =  [11  17  35  31  3  3  5]/100;  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
    c1       =  c0;                  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
    dcr      =  [1,1,1,-1,-1,-1,0]*1e-3;  % amplitude of random noise [wt]
    dcg      =  [-1,-1,-1,1,1,1,0]*1e-2;  % amplitude of centred gaussian [wt]
    dTg      =  5;
    dTr      =  0.01;
    dr_trc   =  [1,1,1,-1,-1,-1].*1e-3;
    dg_trc   =  [-1,-1,-1,1,1,1].*1e-2;

    % closed boundaries for gas flux
    periodic  =  1;
    bndmode   =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
    bnd_w     =  h;                   % boundary layer width [m]
    tau_T     =  bnd_w^2/1e-6;      % wall cooling/assimilation time [s]
    Twall     =  [300,300,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
    Ptop      =  2.0e8;               % top pressure [Pa]
    fin       =  0;
    fout      =  0;

    calID    =  'DEMO';              % phase diagram calibration

    % set numerical model parameters
    TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
    ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
    CFL      =  1;                   % (physical) time stepping courant number (multiplies stable step) [0,1]
    rtol     =  atol/1e6;            % outer its absolute tolerance
    maxit    =  100;                 % maximum outer its
    alpha    =  0.50;                % iterative step size parameter
    Pcouple  =  1;

    % create output directory
    if ~isfolder([outdir,'/',runID])
        mkdir([outdir,'/',runID]);
    end

    % run code
    run('../src/main')

    % plot convergence
    ES = norm(diff(hist.ES(Nt/2:Nt  ))./diff(hist.time(Nt/2:Nt)),'fro')./sqrt(length(diff(hist.ES(Nt/2:Nt))         ));
    EB = norm(diff(hist.EB(Nt/2:Nt  ))./diff(hist.time(Nt/2:Nt)),'fro')./sqrt(length(diff(hist.EB(Nt/2:Nt))         ));
    EM = norm(diff(hist.EM(Nt/2:Nt  ))./diff(hist.time(Nt/2:Nt)),'fro')./sqrt(length(diff(hist.EM(Nt/2:Nt))         ));
    EX = norm(diff(hist.EX(Nt/2:Nt  ))./diff(hist.time(Nt/2:Nt)),'fro')./sqrt(length(diff(hist.EX(Nt/2:Nt))         ));
    EF = norm(diff(hist.EF(Nt/2:Nt  ))./diff(hist.time(Nt/2:Nt)),'fro')./sqrt(length(diff(hist.EF(Nt/2:Nt))         ));
    EC = norm(diff(hist.EC(Nt/2:Nt,:))./repmat(diff(hist.time(Nt/2:Nt)),cal.ncmp,1).','fro')./sqrt(length(diff(hist.EC(Nt/2:Nt))*cal.ncmp));
    ET = norm(diff(hist.ET(Nt/2:Nt,:))./repmat(diff(hist.time(Nt/2:Nt)),cal.ntrc,1).','fro')./sqrt(length(diff(hist.ET(Nt/2:Nt))*cal.ntrc));

    clist = [colororder;[0 0 0]];

    fh20 = figure(20);
    subplot(4,1,1);
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.ES(Nt/2:Nt)))./diff(hist.time(Nt/2:Nt)),'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,2);
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.EB(Nt/2:Nt)))./diff(hist.time(Nt/2:Nt)),'-',LW{:}); hold on; axis tight; box on; hold on
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.EM(Nt/2:Nt)))./diff(hist.time(Nt/2:Nt)),'-',LW{:}); hold on; axis tight; box on;
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.EX(Nt/2:Nt)))./diff(hist.time(Nt/2:Nt)),'-',LW{:}); hold on; axis tight; box on;
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.EF(Nt/2:Nt)))./diff(hist.time(Nt/2:Nt)),'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_{\bar{\rho}}, \Delta E_{F^i}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,3);
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.EC(Nt/2:Nt,:)))./repmat(diff(hist.time(Nt/2:Nt)),cal.ncmp,1).','-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_{C_j}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,4);
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.ET(Nt/2:Nt,:)))./repmat(diff(hist.time(Nt/2:Nt)),cal.ntrc,1).','-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_{\Theta_k}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    xlabel('Time [s]',TX{:},FS{:});

    fh21 = figure(21);
    p1 = loglog(atol,ES,'+','Color',clist(1,:),'MarkerSize',10,'LineWidth',2); hold on; box on;
    p2 = loglog(atol,EB,'s','Color',clist(2,:),'MarkerSize',10,'LineWidth',2);
    p3 = loglog(atol,EM,'o','Color',clist(3,:),'MarkerSize',10,'LineWidth',2);
    p4 = loglog(atol,EX,'d','Color',clist(4,:),'MarkerSize',10,'LineWidth',2);
    p5 = loglog(atol,EF,'*','Color',clist(5,:),'MarkerSize',10,'LineWidth',2);
    p6 = loglog(atol,EC,'^','Color',clist(6,:),'MarkerSize',10,'LineWidth',2);
    p7 = loglog(atol,ET,'v','Color',clist(7,:),'MarkerSize',10,'LineWidth',2);
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('abs. residual tolerance [1]','Interpreter','latex','FontSize',16)
    ylabel('rel. conservation error rate [1/s]','Interpreter','latex','FontSize',16)
    title('Global conservation with nonlinear convergence','Interpreter','latex','FontSize',20)

    if atol == ATOL(end)
        p8 = loglog(ATOL,eps.*ones(size(ATOL)),'k-' ,'LineWidth',2);  % plot trend for comparison
        % legend([p1,p2,p3,p4,p5,p6,p7,p8],{'error $S$','error $\bar{\rho}$','error $C_j$','error $\Theta_k$','machine prec.'},'Interpreter','latex','box','on','location','southeast')
        legend([p1,p2,p3,p4,p5,p6,p7,p8],{'error $S$','error $\bar{\rho}$','error $M$','error $X$','error $F$','error $C_j$','error $\Theta_k$','machine prec.'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

end

name = [outdir,'/',runID,'/',runID,'_',TINT,'_',ADVN];
print(fh21,name,'-dpng','-r300','-vector');