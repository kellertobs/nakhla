clear; close all;

NN = [40,80,160];  % test increasing time steps

for nn = NN
    
% set run parameters
runID    =  'bnchm_mms';         % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  1;                   % output frame plotted/saved every 'nop' time steps
plot_op  =  0;                   % switch on to live plot of results
save_op  =  0;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence
react    =  1;                   % switch on reactive mode
diseq    =  1;                   % switch on disequilibrium approac
bnchm    =  1;                   % switch on to run manufactured solution benchmark on flui mechanics solver

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  nn + 2;              % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
M        =  1;                   % number of time steps to take
hr       =  3600;                % conversion seconds to hours
yr       =  24*365.25*hr;        % conversion seconds to years
tend     =  1*yr;                % end time for simulation [s]
dt       =  10;                  % initial time step [s]
dtmax    =  10;                  % maximum time step [s]

% set initial thermo-chemical state
seed     =  15;                  % random perturbation seed
smth     =  (N/30)^2;            % regularisation of initial random perturbation
zlay     =  0.5;                 % layer thickness (relative to domain depth D)
wlay_T   =  1e-6;                % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0       =  1050;                % temperature top layer [deg C]
T1       =  1050;                % temperature base layer [deg C]
dT       =  0;                   % amplitude of random noise [deg C]
c0       =  0.50;                % major component top layer [wt SiO2]
c1       =  0.50;                % major component base layer [wt SiO2]
dc       =  0e-5;                % amplitude of random noise [wt SiO2]
v0       =  0.03;                % volatile component top layer [wt H2O]
v1       =  0.03;                % volatile component base layer [wt H2O]
dv       =  0e-6;                % amplitude of random noise [wt H2O]

% set model trace and isotope geochemistry parameters
it0      =  1;                   % incompatible tracer top layer [wt ppm]
it1      =  1;                   % incompatible tracer base layer [wt ppm]
dit      =  0.0;                 % incompatible tracer random noise [wt ppm]
KIT      =  1e-2;                % incompatible tracer partition coefficient
ct0      =  1;                   % compatible tracer top layer [wt ppm]
ct1      =  1;                   % compatible tracer base layer [wt ppm]
dct      =  -0.0;                % compatible tracer random noise [wt ppm]
KCT      =  1e2;                 % compatible tracer partition coefficient
si0      =  0;                   % stable isotope ratio top layer [delta]
si1      =  0;                   % stable isotope ratio base layer [delta]
dsi      =  0.0;                 % stable isotope ratio random noise [delta]
ri0      =  1;                   % radiogenic isotope top layer [wt ppm]
ri1      =  1;                   % radiogenic isotope base layer [wt ppm]
dri      = -0.0;                 % radiogenic isotope random noise [wt ppm]
KRIP     =  10;                  % radiogenic parent isotope partition coefficient
KRID     =  0.1;                 % radiogenic daughter isotope partition coefficient
HLRIP    =  1e4*yr;              % radiogenic parent isotope half-life [s]
HLRID    =  1e3*yr;              % radiogenic daughter isotope half-life [s]

% set thermo-chemical boundary parameters
Ptop     =  1e8;                 % top pressure [Pa]
bndmode  =  0;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls)
bndinit  =  0;                   % switch on (1) to initialise with already established boundary layers
dw       =  2*h;                 % boundary layer thickness for cooling/outgassing/assimilation [m]
fin      =  0;                   % ingassing factor (0 = no ingassing; 1 = free flow ingassing)
fout     =  0;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)
tau_T    =  4*hr;                % wall cooling/assimilation time [s]
tau_a    =  8*hr;                % wall cooling/assimilation time [s]
Twall    =  500;                 % wall temperature [degC] (nan = insulating)
cwall    =  nan;                 % wall major component [wt SiO2] (nan = no assimilation)
vwall    =  nan;                 % wall volatile component [wt H2O] (nan = no assimilation)
itwall   =  nan;                 % wall incomp. tracer [wt ppm] (nan = no assimilation)
ctwall   =  nan;                 % wall comp. tracer [wt ppm] (nan = no assimilation)
siwall   =  nan;                 % wall stable isotope [delta] (nan = no assimilation)
riwall   =  nan;                 % wall radiogenic isotope [wt ppm] (nan = no assimilation)

% set thermo-chemical material parameters
kc       =  1e-7;                % chemical diffusivity [m^2/s]
kTm      =  4;                   % melt thermal conductivity [W/m/k]
kTx      =  1;                   % xtal thermal conductivity [W/m/k]
kTf      =  0.02;                % mvp  thermal conductivity [W/m/k]
Cpm      =  1400;                % melt heat capacity [J/kg/K]
Cpx      =  1000;                % xtal heat capacity [J/kg/K]
Cpf      =  2000;                % mvp  heat capacity [J/kg/K]

% set model phase equilibrium parameters
cphs0    =  0.36;                % phase diagram lower bound composition [wt SiO2]
cphs1    =  0.72;                % phase diagram upper bound composition [wt SiO2]
Tphs0    =  750;                 % phase diagram lower bound temperature [degC]
Tphs1    =  1750;                % phase diagram upper bound temperature [degC]
PhDg     =  5.0;                 % Phase diagram curvature factor (> 1)
perCm    =  0.51;                % peritectic liquidus composition [wt SiO2]
perCx    =  0.48;                % peritectic solidus  composition [wt SiO2]
perT     =  1100;                % peritectic temperature [degC]
clap     =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
dTH2O    =  [1300,1000,300];     % solidus shift from water content [degC/wt^0.75]
tau_r    =  60;                  % reaction time [s]
Dsx      = -300;                 % entropy change of crystallisation [J/kg/K]
Dsf      =  400;                 % entropy change of exsolution [J/kg/K]

% set model rheology parameters
etam0    =  100;                 % melt viscosity [Pas]
etaf0    =  0.1;                 % fluid viscosity [Pas]
etax0    =  1e15;                % crystal viscosity [Pas]
Fmc      =  1e+4;                % major component weakening factor of melt viscosity [1]
Fmv      =  0.5;                 % volatile component weakening factor of melt viscosity [1]
Em       =  150e3;               % activation energy melt viscosity [J/mol]
AA       = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
BB       = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
CC       = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% set model buoyancy parameters
rhom0    =  2900;                % melt phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhox0    =  3300;                % crystal phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhof0    =  500;                 % bubble phase ref. density [kg/m3] (at T0,cphs0,Ptop)
aTm      =  3e-5;                % melt thermal expansivity [1/K]
aTx      =  1e-5;                % xtal thermal expansivity [1/K]
aTf      =  1e-4;                % mvp  thermal expansivity [1/K]
gCm      =  0.5;                 % melt compositional expansion [1/wt]
gCx      =  0.6;                 % xtal compositional expansion [1/wt]
bPf      =  1e-8;                % mvp compressibility [1/Pa]
dx       =  1e-9;                % crystal size [m]
df       =  1e-9;                % bubble size [m]
dm       =  1e-9;                % melt film size [m]
g0       =  10.;                 % gravity [m/s2]

% set numerical model parameters
CFL      =  0.25;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'FRM';               % advection scheme ('UPW2', 'UPW3', or 'FRM')
rtol     =  1e-9;                % outer its relative tolerance
atol     =  1e-6;                % outer its absolute tolerance
maxit    =  100;                 % maximum outer its
alpha    =  0.75;                % iterative lag parameter equilibration
beta     =  0.50;                % iterative lag parameter phase diagram
delta    =  20;                  % smoothness of segregation speed
etamin   =  1e1;                 % minimum viscosity for stabilisation
etamax   =  1e7;                 % maximum viscosity for stabilisation
TINY     =  1e-16;               % minimum cutoff phase, component fractions

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
run('../src/main')

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
p1 = loglog(h,errW,'r+','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on; 
p2 = loglog(h,errU,'g+','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on; 
p3 = loglog(h,errP,'b+','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on; 
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
p5 = loglog(dofs,solvetime,'r+','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on; 
set(gca,'TicklabelInterpreter','latex','FontSize',12)
xlabel('\# dofs [1]','Interpreter','latex','FontSize',16)
ylabel('time to solution [s]','Interpreter','latex','FontSize',16)
title('Scaling of direct solver','Interpreter','latex','FontSize',20)

if nn == NN(1)
    p6 = loglog(DOFS,0.95*solvetime*(DOFS./DOFS(1)).^1,'k-','LineWidth',2);  % plot linear trend for comparison
end
if nn == NN(end)
    legend([p5,p6],{'time to solution','linear'},'Interpreter','latex','box','on','location','southeast')
end

end

name = [opdir,'/',runID,'/',runID,'_bnchm'];
print(fh18,name,'-dpng','-r300','-opengl');

name = [opdir,'/',runID,'/',runID,'_sclng'];
print(fh19,name,'-dpng','-r300','-opengl');
