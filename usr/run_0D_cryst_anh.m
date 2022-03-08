clear; close all;

% set run parameters
runID    =  '0D_cryst_anh';      % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  1;                   % switch on to live plot iterative convergence
react    =  1;                   % switch on reactive mode
diseq    =  1;                   % switch on disequilibrium approac
bnchm    =  0;                   % switch on to run manufactured solution benchmark on flui mechanics solver

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  1 + 2;               % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
M        =  5e4;                 % number of time steps to take
hr       =  3600;                % conversion seconds to hours
yr       =  24*365.25*hr;        % conversion seconds to years
tend     =  12*hr;               % end time for simulation [s]
dt       =  5;                   % initial time step [s]
dtmax    =  5;                   % maximum time step [s]

% set initial thermo-chemical state
seed     =  15;                  % random perturbation seed
smth     =  (N/30)^2;            % regularisation of initial random perturbation
zlay     =  0.5;                 % layer thickness (relative to domain depth D)
wlay_T   =  1e-6;                % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0       =  1350;                % temperature top layer [deg C]
T1       =  1350;                % temperature base layer [deg C]
dT       =  0;                   % amplitude of random noise [deg C]
c0       =  0.49;                % major component top layer [wt SiO2]
c1       =  0.49;                % major component base layer [wt SiO2]
dc       =  0e-5;                % amplitude of random noise [wt SiO2]
v0       =  0.00;                % volatile component top layer [wt H2O]
v1       =  0.00;                % volatile component base layer [wt H2O]
dv       =  0e-6;                % amplitude of random noise [wt H2O]

% set model trace and isotope geochemistry parameters
it0      =  1;                   % incompatible tracer top layer [wt ppm]
it1      =  1;                   % incompatible tracer base layer [wt ppm]
dit      =  0.00;                % incompatible tracer random noise [wt ppm]
KIT      =  1e-2;                % incompatible tracer partition coefficient
ct0      =  1;                   % compatible tracer top layer [wt ppm]
ct1      =  1;                   % compatible tracer base layer [wt ppm]
dct      =  -0.00;               % compatible tracer random noise [wt ppm]
KCT      =  1e2;                 % compatible tracer partition coefficient
si0      =  -1;                  % stable isotope ratio top layer [delta]
si1      =  1;                   % stable isotope ratio base layer [delta]
dsi      =  0.00;                % stable isotope ratio random noise [delta]
ri0      =  1;                   % radiogenic isotope top layer [wt ppm]
ri1      =  1;                   % radiogenic isotope base layer [wt ppm]
dri      =  -0.00;               % radiogenic isotope random noise [wt ppm]
KRIP     =  10;                  % radiogenic parent isotope partition coefficient
KRID     =  0.1;                 % radiogenic daughter isotope partition coefficient
HLRIP    =  1e4*yr;              % radiogenic parent isotope half-life [s]
HLRID    =  1e3*yr;              % radiogenic daughter isotope half-life [s]

% set thermo-chemical boundary parameters
Ptop     =  1e8;                 % top pressure [Pa]
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls)
bndinit  =  0;                   % switch on (1) to initialise with already established boundary layers
dw       =  2*h;                 % boundary layer thickness for assimilation [m]
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
tau_r    =  10;                  % reaction time [s]
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
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]
dm       =  1e-4;                % melt film size [m]
g0       =  10.;                 % gravity [m/s2]

% set numerical model parameters
CFL      =  0.5;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'FRM';               % advection scheme ('UPW2', 'UPW3', or 'FRM')
rtol     =  1e-5;                % outer its relative tolerance
atol     =  1e-7;                % outer its absolute tolerance
maxit    =  100;                 % maximum outer its
alpha    =  0.75;                % iterative lag parameter equilibration
beta     =  0.50;                % iterative lag parameter phase diagram
delta    =  20;                  % smoothness of segregation speed
etamin   =  1e1;                 % minimum viscosity for stabilisation
etamax   =  1e8;                 % maximum viscosity for stabilisation
TINY     =  1e-16;               % minimum cutoff phase, component fractions

% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 
    parfile = [opdir,'/',runID,'/',runID,'_par'];
    save(parfile);
end

% run code
run('../src/main')
