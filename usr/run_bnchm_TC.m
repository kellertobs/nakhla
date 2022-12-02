clear; close all;

DT = [1/4,1/8,1/16];

for dt = DT
    
% set run parameters
runID    =  'bnchm_cnsv';        % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  4;                   % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
save_op  =  0;                   % switch on to save output to file
plot_cv  =  1;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  100 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
M        =  4/dt*DT(1);          % number of time steps to take
hr       =  3600;                % conversion seconds to hours
yr       =  24*365.25*hr;        % conversion seconds to years
tend     =  1*yr;                % end time for simulation [s]
dtmax    =  dt;                  % maximum time step [s]

% set initial thermo-chemical state
seed     =  15;                  % random perturbation seed
smth     =  10;                  % regularisation of initial random perturbation
zlay     =  0.5;                 % layer thickness (relative to domain depth D)
wlay_T   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0       =  1055;                % temperature top layer [deg C]
T1       =  1055;                % temperature base layer [deg C]
dT       =  0;                   % amplitude of random noise [deg C]
c0       =  0.52;                % major component top layer [wt SiO2]
c1       =  0.52;                % major component base layer [wt SiO2]
dc       =  1e-3;                % amplitude of random noise [wt SiO2]
v0       =  0.04;                % volatile component top layer [wt H2O]
v1       =  0.04;                % volatile component base layer [wt H2O]
dv       =  1e-4;                % amplitude of random noise [wt H2O]

% set model trace and isotope geochemistry parameters
te0      =  [1,1,1,1];           % trace elements top layer [wt ppm]
te1      =  [1,1,1,1];           % trace elements base layer [wt ppm]
dte      =  1e-3.*[-1,-1,1,1];   % trace elements random noise [wt ppm]
Kte      =  [0.01,0.1,3,10];     % trace elements partition coefficients
ir0      =  [0,-1];              % isotope ratios top layer [delta]
ir1      =  [0, 1];              % isotope ratios base layer [delta]
dir      =  [1, 0];              % isotope ratios random noise [delta]

% set thermo-chemical boundary parameters
Ptop     =  1.25e8;              % top pressure [Pa]
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls)
bndinit  =  0;                   % switch on (1) to initialise with internal boundary layers
dw       =  1*h;                 % boundary layer thickness [m]
fin      =  0;                   % ingassing factor (0 = no ingassing; 1 = free flow ingassing)
fout     =  0;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)
tau_T    =  8*hr;                % wall cooling/assimilation time [s]
tau_a    =  4*hr;                % wall cooling/assimilation tie [s]
Twall    =  300;                 % wall temperature [degC] (nan = insulating)
cwall    =  nan;                 % wall major component [wt SiO2] (nan = no assimilation)
vwall    =  nan;                 % wall volatile component [wt H2O] (nan = no assimilation)
tewall   =  [nan,nan,nan,nan];   % wall trace elements [wt ppm] (nan = no assimilation)
irwall   =  [nan,nan,nan,nan];   % wall isotope ratios [delta] (nan = no assimilation)

% set thermo-chemical material parameters
calID    =  'morb';              % phase diagram calibration
kT0      =  4;                   % thermal conductivity [W/m/K]
cP       =  1200;                % heat capacity [J/kg/K]
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  400;                 % entropy change of exsolution [J/kg]

% set model buoyancy parameters
d0       =  1e-3;                % crystal/bubble size [m]
g0       =  10.;                 % gravity [m/s2]

% set numerical model parameters
CFL      =  1.00;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-9;                % outer its absolute tolerance
maxit    =  50;                  % maximum outer its
lambda   =  0.50;                % iterative lag parameter equilibration
etareg   =  1e0;                 % viscosity regularisation parameter

% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

% save input parameters and runtime options (unless restarting or benchmarking)
% if restart == 0 && bnchm == 0
%     parfile = [opdir,'/',runID,'/',runID,'_par'];
%     save(parfile);
% end

% run code
run('../src/main')

run('../src/output')

% plot convergence
EM = norm(diff(hist.EM(2:end)),'fro')./sqrt(length(diff(hist.EM(2:end))));
ES = norm(diff(hist.ES(2:end)),'fro')./sqrt(length(diff(hist.ES(2:end))));
EC = norm(diff(hist.EC(2:end)),'fro')./sqrt(length(diff(hist.EC(2:end))));
EV = norm(diff(hist.EV(2:end)),'fro')./sqrt(length(diff(hist.EV(2:end))));
EO = norm(diff(hist.EC_oxd(2:end,:)),'fro')./sqrt(length(reshape(diff(hist.EC_oxd(2:end,:)),1,[])));

fh15 = figure(15);
p1 = loglog(dt,EM,'kd','MarkerSize',8,'LineWidth',2); hold on; box on;
p2 = loglog(dt,ES,'rs','MarkerSize',8,'LineWidth',2);
p3 = loglog(dt,EC,'go','MarkerSize',8,'LineWidth',2);
p4 = loglog(dt,EV,'bv','MarkerSize',8,'LineWidth',2);
p5 = loglog(dt,EO,'m+','MarkerSize',8,'LineWidth',2);
set(gca,'TicklabelInterpreter','latex','FontSize',12)
xlabel('time step [s]','Interpreter','latex','FontSize',16)
ylabel('rel. numerical error [1]','Interpreter','latex','FontSize',16)
title('Numerical convergence in time','Interpreter','latex','FontSize',20)

if dt == DT(1)
    p6 = loglog(DT,geomean([EC,EV,ES]).*(DT./DT(1)).^1,'k-' ,'LineWidth',2);  % plot trend for comparison
    p6 = loglog(DT,geomean([EC,EV,ES]).*(DT./DT(1)).^2,'k-' ,'LineWidth',2);  % plot trend for comparison
    p6 = loglog(DT,geomean([EC,EV,ES]).*(DT./DT(1)).^3,'k-' ,'LineWidth',2);  % plot trend for comparison
end
if dt == DT(end)
    legend([p1,p2,p3,p4,p5,p6],{'error $M$','error $S$','error $C$','error $V$','error $C_{oxd}$','trends'},'Interpreter','latex','box','on','location','southeast')
end
drawnow;

end

name = [opdir,'/',runID,'/',runID,'_bnchm'];
print(fh15,name,'-dpng','-r300','-vector');