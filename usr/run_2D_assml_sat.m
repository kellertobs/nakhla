clear all; close all;

% set run parameters
runID    =  '2D_assml_sat';      % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence
react    =  1;                   % switch on reactive mode
diseq    =  1;                   % switch on disequilibrium approac
bnchm    =  0;                   % switch on to run manufactured solution benchmark on flui mechanics solver

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  200 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
M        =  1e5;                 % number of time steps to take
hr       =  3600;                % conversion seconds to hours
yr       =  24*365.25*hr;        % conversion seconds to years
tend     =  1*yr;                % end time for simulation [s]
dt       =  100;                 % initial time step [s]
dtmax    =  1e3;                 % maximum time step [s]

% set initial thermo-chemical state
seed     =  15;                  % random perturbation seed
smth     =  (N/30)^2;            % regularisation of initial random perturbation
zlay     =  0.5;                 % layer thickness (relative to domain depth D)
wlay_T   =  4*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0       =  1250;                % temperature top layer [deg C]
T1       =  1250;                % temperature base layer [deg C]
dT       =  0;                   % amplitude of random noise [deg C]
c0       =  0.49;                % major component top layer [wt SiO2]
c1       =  0.49;                % major component base layer [wt SiO2]
dc       =  1e-4;                % amplitude of random noise [wt SiO2]
v0       =  0.04;                % volatile component top layer [wt H2O]
v1       =  0.04;                % volatile component base layer [wt H2O]
dv       =  1e-4;                % amplitude of random noise [wt H2O]

% set model trace and isotope geochemistry parameters
it0      =  0.1;                 % incompatible tracer top layer [wt ppm]
it1      =  0.1;                 % incompatible tracer base layer [wt ppm]
dit      =  0.00;                % incompatible tracer random noise [wt ppm]
KIT      =  1e-2;                % incompatible tracer partition coefficient
ct0      =  10.;                 % compatible tracer top layer [wt ppm]
ct1      =  10.;                 % compatible tracer base layer [wt ppm]
dct      = -0.00;                % compatible tracer random noise [wt ppm]
KCT      =  1e+2;                % compatible tracer partition coefficient
si0      =  5.;                  % stable isotope ratio top layer [delta]
si1      =  5.;                  % stable isotope ratio base layer [delta]
dsi      =  0.00;                % stable isotope ratio random noise [delta]
ri0      =  2.0;                 % radiogenic isotope top layer [wt ppm]
ri1      =  2.0;                 % radiogenic isotope base layer [wt ppm]
dri      = -0.00;                % radiogenic isotope random noise [wt ppm]
KRIP     =  10.;                 % radiogenic parent isotope partition coefficient
KRID     =  0.1;                 % radiogenic daughter isotope partition coefficient
HLRIP    =  100*yr;              % radiogenic parent isotope half-life [s]
HLRID    =    1*yr;              % radiogenic daughter isotope half-life [s]

% set thermo-chemical boundary parameters
Ptop     =  1e8;                 % top pressure [Pa]
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls)
bndinit  =  0;                   % switch on (1) to initialise with already established boundary layers
dw       =  1*h;                 % boundary layer thickness for assimilation [m]
fin      =  1;                   % ingassing factor (0 = no ingassing; 1 = free flow ingassing)
fout     =  1;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)
tau_T    =  4*hr;                % wall cooling/assimilation time [s]
tau_a    =  4*hr;                % wall cooling/assimilation time [s]
Twall    =  500;                 % wall temperature [degC] (nan = insulating)
cwall    =  0.70;                % wall major component [wt SiO2] (nan = no assimilation)
vwall    =  0.02;                % wall volatile component [wt H2O] (nan = no assimilation)
itwall   =  10;                  % wall incomp. tracer [wt ppm] (nan = no assimilation)
ctwall   =  0.1;                 % wall comp. tracer [wt ppm] (nan = no assimilation)
siwall   =  -5;                  % wall stable isotope [delta] (nan = no assimilation)
riwall   =  0.5;                 % wall radiogenic isotope [wt ppm] (nan = no assimilation)

% set thermo-chemical material parameters
kc       =  1e-4;                % chemical diffusivity [kg/m/s]
kTm      =  4;                   % melt thermal conductivity [W/m/K]
kTx      =  1;                   % xtal thermal conductivity [W/m/K]
kTf      =  0.02;                % mvp  thermal conductivity [W/m/K]
Cp       =  1300;                % heat capacity [J/kg/K]
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  500;                 % entropy change of exsolution [J/kg]

% set phase diagram parameters
cphs0    =  0.43;                % phase diagram lower bound composition [wt SiO2]
cphs1    =  0.78;                % phase diagram upper bound composition [wt SiO2]
Tphs0    =  844;                 % phase diagram lower bound temperature [degC]
Tphs1    =  1867;                % phase diagram upper bound temperature [degC]
PhDg     =  [7.0,4.2,1.0,0.9];   % Phase diagram curvature factor (> 1)
perCm    =  0.517;               % peritectic liquidus composition [wt SiO2]
perCx    =  0.475;               % peritectic solidus  composition [wt SiO2]
perT     =  1147;                % peritectic temperature [degC]
clap     =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
dTH2O    =  [1200,1000,100];     % solidus shift from water content [degC/wt^0.75]
tau_r    =  60;                  % reaction time [s]

% set model rheology parameters
etam0    =  300;                 % melt viscosity [Pas]
etaf0    =  0.1;                 % fluid viscosity [Pas]
etax0    =  1e15;                % crystal viscosity [Pas]
Fmc      =  1e+5;                % major component weakening factor of melt viscosity [1]
Fmv      =  0.6;                 % volatile component weakening factor of melt viscosity [1]
Em       =  175e3;               % activation energy melt viscosity [J/mol]
AA       = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
BB       = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
CC       = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% set model buoyancy parameters
rhom0    =  2750;                % melt phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhox0    =  3050;                % crystal phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhof0    =  1000;                % bubble phase ref. density [kg/m3] (at T0,cphs0,Ptop)
aTm      =  3e-5;                % melt thermal expansivity [1/K]
aTx      =  1e-5;                % xtal thermal expansivity [1/K]
aTf      =  1e-4;                % mvp  thermal expansivity [1/K]
gCm      =  0.5;                 % melt compositional expansion [1/wt]
gCx      =  0.5;                 % xtal compositional expansion [1/wt]
bPf      =  1e-8;                % mvp compressibility [1/Pa]
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]
g0       =  10.;                 % gravity [m/s2]

% set numerical model parameters
CFL      =  0.25;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'FRM';               % advection scheme ('UPW2', 'UPW3', or 'FRM')
theta    =  0.5;                 % time-stepping parameter (1 = 1st-order implicit; 1/2 = 2nd-order semi-implicit)
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-7;                % outer its absolute tolerance
maxit    =  20;                  % maximum outer its
alpha    =  0.25;                % iterative lag parameter equilibration
beta     =  0.10;                % iterative lag parameter phase diagram
delta    =  2;                   % smoothness of segregation speed
etamin   =  1e2;                 % minimum viscosity for stabilisation
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
