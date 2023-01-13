% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% test decreasing time step
RTOL = [1e-3,1e-6,1e-9];

for rtol = RTOL

    % set run parameters
    runID    =  'bnchm_cnsv';        % run identifier
    nop      =  10;                  % output frame plotted/saved every 'nop' time steps
    plot_op  =  1;                   % switch on to live plot of results
    plot_cv  =  1;                   % switch on to live plot iterative convergence

    % set model timing parameters
    Nt       =  nop;                 % number of time steps to take
    dt       =  1/4;                 % set initial time step

    % set initial thermo-chemical state
    T0       =  1140;                % temperature top layer [deg C]
    c0       =  0.51;                % major component top layer [wt SiO2]
    dcr      =  1e-4;                % amplitude of random noise [wt SiO2]
    v0       =  0.04;                % volatile component top layer [wt H2O]
    dvr      =  1e-5;                % amplitude of random noise [wt H2O]

    % closed boundaries for gas flux
    fin      =  0;
    fout     =  0;

    % set numerical model parameters
    CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
    atol     =  rtol/1e3;            % outer its absolute tolerance
    maxit    =  50;                  % maximum outer its
    lambda   =  0.50;                % iterative lag parameter equilibration

    % create output directory
    if ~isfolder([opdir,'/',runID])
        mkdir([opdir,'/',runID]);
    end

    % run code
    run('../src/main')

    % plot convergence
    EM = norm(diff(hist.EM(Nt/2:Nt))/dt,'fro')./sqrt(length(diff(hist.EM(Nt/2:Nt))));
    ES = norm(diff(hist.ES(Nt/2:Nt))/dt,'fro')./sqrt(length(diff(hist.ES(Nt/2:Nt))));
    EC = norm(diff(hist.EC(Nt/2:Nt))/dt,'fro')./sqrt(length(diff(hist.EC(Nt/2:Nt))));
    EV = norm(diff(hist.EV(Nt/2:Nt))/dt,'fro')./sqrt(length(diff(hist.EV(Nt/2:Nt))));

    fh14 = figure(14);
    subplot(4,1,1);
    plot((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,diff(hist.EM(Nt/2:Nt))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,2);
    plot((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,diff(hist.ES(Nt/2:Nt))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,3);
    plot((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,diff(hist.EC(Nt/2:Nt))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_C$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,4);
    plot((Nt/2+1/2)*dt:dt:(Nt-1/2)*dt,diff(hist.EV(Nt/2:Nt))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_V$',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel('Time [s]',TX{:},FS{:});

    fh15 = figure(15);
    p1 = loglog(rtol,EM,'kd','MarkerSize',8,'LineWidth',2); hold on; box on;
    p2 = loglog(rtol,ES,'rs','MarkerSize',8,'LineWidth',2);
    p3 = loglog(rtol,EC,'go','MarkerSize',8,'LineWidth',2);
    p4 = loglog(rtol,EV,'bv','MarkerSize',8,'LineWidth',2);
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('rel. nonlinear tolerance [1]','Interpreter','latex','FontSize',16)
    ylabel('rel. numerical error rate [1/s]','Interpreter','latex','FontSize',16)
    title('Global conservation in time','Interpreter','latex','FontSize',20)

    if rtol == RTOL(1)
        p5 = loglog(RTOL,1e-16.*ones(size(RTOL)),'k-' ,'LineWidth',2);  % plot trend for comparison
    end
    if rtol == RTOL(end)
        legend([p1,p2,p3,p4,p5],{'error $M$','error $S$','error $C$','error $V$','machine prec.'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

end

name = [opdir,'/',runID,'/',runID,'_bnchm'];
print(fh15,name,'-dpng','-r300','-vector');