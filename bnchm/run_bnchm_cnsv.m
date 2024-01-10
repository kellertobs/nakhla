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
    N        =  100;                 % number of grid points in z-direction (incl. 2 ghosts)
    h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

    % set model timing parameters
    Nt       =  nop;                 % number of time steps to take
    dt       =  10;                  % set initial time step

    % set initial thermo-chemical state
    T0       =  1050;                % temperature top  layer [deg C]
    T1       =  T0;                  % temperature base layer [deg C]
    c0       =  [0.04  0.12  0.44  0.24  0.14  0.02  0.04];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
    c1       =  c0;                  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
    dcr      =  [1,1,1,-1,-1,-1,0]*1e-4;  % amplitude of random noise [wt]
    dcg      =  [-1,-1,-1,1,1,1,0]*1e-3;  % amplitude of centred gaussian [wt]

    % closed boundaries for gas flux
    fin      =  0;
    fout     =  0;

    % set numerical model parameters
    TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
    ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
    CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
    atol     =  rtol/1e6;            % outer its absolute tolerance
    maxit    =  100;                 % maximum outer its
    alpha    =  0.75;                % iterative step size parameter
    beta     =  0.10;                % iterative damping parameter

    % create output directory
    if ~isfolder([opdir,'/',runID])
        mkdir([opdir,'/',runID]);
    end

    % run code
    run('../src/main')

    % plot convergence
    EM = norm(diff(hist.EM(Nt/2:Nt  ))/dt,'fro')./sqrt(length(diff(hist.EM(Nt/2:Nt))         ));
    ES = norm(diff(hist.ES(Nt/2:Nt  ))/dt,'fro')./sqrt(length(diff(hist.ES(Nt/2:Nt))         ));
    EC = norm(diff(hist.EC(Nt/2:Nt,:))/dt,'fro')./sqrt(length(diff(hist.EC(Nt/2:Nt))*cal.ncmp));

    fh14 = figure(14);
    subplot(3,1,1);
    semilogy((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,diff(hist.EM(Nt/2:Nt))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(3,1,2);
    semilogy((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,diff(hist.ES(Nt/2:Nt))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(3,1,3);
    semilogy((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,diff(hist.EC(Nt/2:Nt,:))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_C$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    xlabel('Time [s]',TX{:},FS{:});

    fh15 = figure(15);
    p1 = loglog(rtol,EM,'kd','MarkerSize',8,'LineWidth',2); hold on; box on;
    p2 = loglog(rtol,ES,'rs','MarkerSize',8,'LineWidth',2);
    p3 = loglog(rtol,EC,'go','MarkerSize',8,'LineWidth',2);
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('rel. nonlinear tolerance [1]','Interpreter','latex','FontSize',16)
    ylabel('rel. numerical error rate [1/s]','Interpreter','latex','FontSize',16)
    title('Global conservation in time','Interpreter','latex','FontSize',20)

    if rtol == RTOL(end)
        p5 = loglog(RTOL,1e-16.*ones(size(RTOL)),'k-' ,'LineWidth',2);  % plot trend for comparison
        legend([p1,p2,p3,p5],{'error $M$','error $S$','error $C$','machine prec.'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

end

name = [opdir,'/',runID,'/',runID,'_',TINT];
print(fh15,name,'-dpng','-r300','-vector');