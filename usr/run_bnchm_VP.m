clear; close all;

NN = [50,100,200];  % test increasing time steps

for nn = NN
    
% set run parameters
runID    =  'bnchm_mms';         % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  1;                   % output frame plotted/saved every 'nop' time steps
plot_op  =  0;                   % switch on to live plot of results
save_op  =  0;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence
bnchm    =  1;                   % switch on to run manufactured solution benchmark on flui mechanics solver

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  nn + 2;              % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
M        =  1e5;                 % number of time steps to take
hr       =  3600;                % conversion seconds to hours
yr       =  24*365.25*hr;        % conversion seconds to years
tend     =  1*yr;                % end time for simulation [s]
dt       =  0;                   % initial time step [s]
dtmax    =  0;                   % maximum time step [s]

% set initial thermo-chemical state
seed     =  15;                  % random perturbation seed
smth     =  (N/30)^2;            % regularisation of initial random perturbation
zlay     =  0.5;                 % layer thickness (relative to domain depth D)
wlay_T   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0       =  1100;                % temperature top layer [deg C]
T1       =  1100;                % temperature base layer [deg C]
dT       =  0;                   % amplitude of random noise [deg C]
c0       =  0.52;                % major component top layer [wt SiO2]
c1       =  0.52;                % major component base layer [wt SiO2]
dc       =  1e-5;                % amplitude of random noise [wt SiO2]
v0       =  0.04;                % volatile component top layer [wt H2O]
v1       =  0.04;                % volatile component base layer [wt H2O]
dv       =  0e-6;                % amplitude of random noise [wt H2O]

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
fout     =  1;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)
tau_T    =  10*hr;               % wall cooling/assimilation time [s]
tau_a    =  1*hr;                % wall cooling/assimilation tie [s]
Twall    =  300;                 % wall temperature [degC] (nan = insulating)
cwall    =  nan;                 % wall major component [wt SiO2] (nan = no assimilation)
vwall    =  nan;                 % wall volatile component [wt H2O] (nan = no assimilation)
tewall   =  [nan,nan,nan,nan];   % wall trace elements [wt ppm] (nan = no assimilation)
irwall   =  [nan,nan,nan,nan];   % wall isotope ratios [delta] (nan = no assimilation)

% set thermo-chemical material parameters
calID    =  'morb';             % phase diagram calibration
kT0      =  4;                   % thermal conductivity [W/m/K]
cP       =  1200;                % heat capacity [J/kg/K]
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  400;                 % entropy change of exsolution [J/kg]

% set model buoyancy parameters
d0       =  1e-3;                % crystal/bubble size [m]
g0       =  10.;                 % gravity [m/s2]

% set numerical model parameters
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
BCA      =  {'',''};             % boundary condition on advection (top/bot, sides)
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-7;                % outer its absolute tolerance
maxit    =  10;                  % maximum outer its
lambda   =  0.25;                % iterative lag parameter equilibration
etareg   =  1e0;                 % viscosity regularisation parameter

% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

% save input parameters and runtime options (unless restarting or benchmarking)
if restart == 0 && bnchm == 0
    parfile = [opdir,'/',runID,'/',runID,'_par'];
    save(parfile);
end

% run code
run('../src/mms')

figure(17); clf;
colormap(ocean);
subplot(2,3,1); imagesc(x_mms,zw_mms,-W(:,2:end-1)*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,2); imagesc(xu_mms,z_mms, U(2:end-1,:)*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,3); imagesc(x_mms ,z_mms, P(2:end-1,2:end-1)/1e3); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $P$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,4); imagesc(x_mms,zw_mms,-(W(:,2:end-1)-W_mms(:,2:end-1))*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,5); imagesc(xu_mms,z_mms, (U(2:end-1,:)-U_mms(2:end-1,:))*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,6); imagesc(x_mms ,z_mms, (P(2:end-1,2:end-1)-P_mms(2:end-1,2:end-1))/1e3); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $P$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
drawnow;

% get solution error
errW = norm(W(:,2:end-1)-W_mms(:,2:end-1),2)./norm(W_mms(:,2:end-1),2);
errU = norm(U(2:end-1,:)-U_mms(2:end-1,:),2)./norm(U_mms(2:end-1,:),2);
errP = norm(P(2:end-1,2:end-1)-P_mms(2:end-1,2:end-1),2)./norm(P_mms(2:end-1,2:end-1),2);

% plot error convergence
fh18 = figure(18); 
p1 = loglog(h,errW,'rs','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on; 
p2 = loglog(h,errU,'go','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on; 
p3 = loglog(h,errP,'bv','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on; 
set(gca,'TicklabelInterpreter','latex','FontSize',12)
xlabel('grid spacing [m]','Interpreter','latex')
ylabel('numerical error [1]','Interpreter','latex')
set(gca,'TicklabelInterpreter','latex')
title('Numerical convergence in space','Interpreter','latex','FontSize',20)
    
if nn == NN(1)
    p4 = loglog(L./NN,mean([errW,errU,errP]).*(NN(1)./NN).^2,'k-','LineWidth',2);  % plot linear trend for comparison
end
if nn == NN(end)
    legend([p1,p2,p3,p4],{'error W','error U','error P','quadratic'},'Interpreter','latex','box','on','location','southeast')
end

% plot error convergence
fh19 = figure(19);
DOFS = (NN+2).*(NN+2) + 2.*(NN+1).*(NN+2);
dofs = (nn+2).*(nn+2) + 2.*(nn+1).*(nn+2);
p5 = loglog(dofs,FMtime,'r+','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on; 
set(gca,'TicklabelInterpreter','latex','FontSize',12)
xlabel('\# dofs [1]','Interpreter','latex','FontSize',16)
ylabel('time to solution [s]','Interpreter','latex','FontSize',16)
title('Scaling of direct solver','Interpreter','latex','FontSize',20)

if nn == NN(1)
    p6 = loglog(DOFS,0.95*FMtime*(DOFS./DOFS(1)).^1,'k-','LineWidth',2);  % plot linear trend for comparison
    p7 = loglog(DOFS,0.95*FMtime*(DOFS./DOFS(1)).^2,'k--','LineWidth',2);  % plot quadratic trend for comparison
end
if nn == NN(end)
    legend([p5,p6,p7],{'time to solution','linear','quadratic'},'Interpreter','latex','box','on','location','southeast')
end

end

name = [opdir,'/',runID,'/',runID,'_bnchm'];
print(fh18,name,'-dpng','-r300','-vector');

name = [opdir,'/',runID,'/',runID,'_sclng'];
print(fh19,name,'-dpng','-r300','-vector');
