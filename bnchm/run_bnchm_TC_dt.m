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
save_op  =  0;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  100 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set initial thermo-chemical state
T0       =  1140;                % temperature top layer [deg C]
dTg      =  80;                  % amplitude of centred gaussian [deg C]
c0       =  0.51;                % major component top layer [wt SiO2]
dcg      = -0.01;                % amplitude of centred gaussian [wt SiO2]
v0       =  0.04;                % volatile component top layer [wt H2O]
dvg      =  0.005;               % amplitude of centred gaussian [wt H2O]

% set model trace and isotope geochemistry parameters (must match # trace elements and isotope ratios in calibration!)
te0      =  [1,1,1,1];           % trace elements top layer [wt ppm]
dteg     =  [1,1,1,1];           % trace elements centred gaussian [wt ppm]
ir0      =  [1, 1];              % isotope ratios top layer [delta]
dirg     =  [1, 1];              % isotope ratios centred gaussian [delta]

fin = 0; fout = 0; Twall = nan;

% set numerical model parameters
CFL      =  1.00;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-9;                % outer its absolute tolerance
maxit    =  50;                  % maximum outer its

% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

cd ../src

DT = [h/2,h/4,h/8];

for dti = DT
    
    dt    =  dti;
    dtmax =  dti;
    Nt    =  L/dti;

    % initialise fields
    init;

    % set velocities to constant values for lateral translation with no segregation
    W(:) = 0; Wm(:) = 0; Wx(:) = 0; Wf(:) = 0; wx(:) = 0;  wm(:) = 0;  wf(:) = 0;
    U(:) = 0; Um(:) = 1; Ux(:) = 1; Uf(:) = 1;
    P(:) = 0;

    % set periodic boundary conditions for advection
    BCA = {'closed','periodic'};

    % set diffusion parameters to zero to isolate advection
    ks(:) = 0; kc(:) = 0; kv(:) = 0;

    % set parameters for non-dissipative, non-reactive flow
    diss(:) = 0;
    tau_r   = 1e32;
    resnorm_VP = 0;

    rhoin = rho;
    Sin = S;
    Cin = C;
    Vin = V;
    Xin = X;
    Fin = F;

    figure(100); clf;
        plot(XX(N/2,inx),Xin(N/2,inx)./rhoin(N/2,inx),'k',XX(N/2,inx),X(N/2,inx)./rho(N/2,inx),'r','LineWidth',1.5); axis tight; box on;
    drawnow;

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

        if step==1
            alpha1 = 1; alpha2 = 1; alpha3 = 0;
            beta1  = 1; beta2  = 0; beta3  = 0;
        elseif step==2
            alpha1 = 1;   alpha2 = 1;   alpha3 = 0;
            beta1  = 1/2; beta2  = 1/2; beta3  = 0;
        else
            alpha1 = 3/2; alpha2 = 4/2; alpha3 = -1/2;
            beta1  = 3/4; beta2  = 2/4; beta3  = -1/4;
        end

        % store previous solution
        Soo = So; So = S;
        Coo = Co; Co = C;
        Voo = Vo; Vo = V;
        Xoo = Xo; Xo = X;
        Foo = Fo; Fo = F;
        TEoo = TEo; TEo = TE;
        IRoo = IRo; IRo = IR;
        dSdtoo = dSdto; dSdto = dSdt;
        dCdtoo = dCdto; dCdto = dCdt;
        dVdtoo = dVdto; dVdto = dVdt;
        dXdtoo = dXdto; dXdto = dXdt;
        dFdtoo = dFdto; dFdto = dFdt;
        dTEdtoo = dTEdto; dTEdto = dTEdt;
        dIRdtoo = dIRdto; dIRdto = dIRdt;
        rhooo = rhoo; rhoo = rho;
        Div_rhoVoo = Div_rhoVo; Div_rhoVo = Div_rhoV;
        Div_Vo  = Div_V;
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
            diss(:) = 0;

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
        plot(XX(N/2,inx),Xin(N/2,inx)./rhoin(N/2,inx),'k',XX(N/2,inx),X(N/2,inx)./rho(N/2,inx),'r','LineWidth',1.5); axis tight; box on;
        drawnow;
    end


    % plot convergence
    EM = norm(rho-rhoin,'fro')./norm(rhoin,'fro');
    ES = norm(S-Sin,'fro')./norm(Sin,'fro');
    EC = norm(C-Cin,'fro')./norm(Cin,'fro');
    EV = norm(V-Vin,'fro')./norm(Vin,'fro');
    EX = norm(X-Xin,'fro')./norm(Xin,'fro');
    EF = norm(F-Fin,'fro')./norm(Fin,'fro');

    fh15 = figure(15);
    p1 = loglog(dt,EM,'kd','MarkerSize',8,'LineWidth',2); hold on; box on;
    p2 = loglog(dt,ES,'rs','MarkerSize',8,'LineWidth',2);
    p3 = loglog(dt,EC,'go','MarkerSize',8,'LineWidth',2);
    p4 = loglog(dt,EV,'bv','MarkerSize',8,'LineWidth',2);
    p5 = loglog(dt,EX,'m+','MarkerSize',8,'LineWidth',2);
    p6 = loglog(dt,EF,'c^','MarkerSize',8,'LineWidth',2);
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('time step [s]','Interpreter','latex','FontSize',16)
    ylabel('rel. numerical error [1]','Interpreter','latex','FontSize',16)
    title('Numerical convergence in time','Interpreter','latex','FontSize',20)

    if dt == DT(1)
        p7 = loglog(DT,geomean([EM,ES,EC,EV,EX,EF]).*(DT./DT(1)).^1,'k-' ,'LineWidth',2);  % plot trend for comparison
        p8 = loglog(DT,geomean([EM,ES,EC,EV,EX,EF]).*(DT./DT(1)).^2,'k--','LineWidth',2);  % plot trend for comparison
    end
    if dt == DT(end)
        legend([p1,p2,p3,p4,p5,p6,p7,p8],{'error $M$','error $S$','error $C$','error $V$','error $X$','error $F$','linear','quadratic'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

end

name = [opdir,'/',runID,'/',runID,'_bnchm'];
print(fh15,name,'-dpng','-r300','-vector');