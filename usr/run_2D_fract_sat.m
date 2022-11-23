% prep workspace
clear all; close all;

% set run parameters
runID    =  '2D_fract_sat';      % run identifier
opdir    =  '../out';            % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  200;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence
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
dt       =  1;                   % initial time step [s]
dtmax    =  100;                 % maximum time step [s]

% set initial thermo-chemical state
seed     =  15;                  % random perturbation seed
smth     =  (N/30)^2;            % regularisation of initial random perturbation
zlay     =  0.5;                 % layer thickness (relative to domain depth D)
wlay_T   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0       =  1050;                % temperature top layer [deg C]
T1       =  1050;                % temperature base layer [deg C]
dT       =  0;                   % amplitude of random noise [deg C]
c0       =  0.52;                % major component top layer [wt SiO2]
c1       =  0.52;                % major component base layer [wt SiO2]
dc       =  1e-5;                % amplitude of random noise [wt SiO2]
v0       =  0.04;                % volatile component top layer [wt H2O]
v1       =  0.04;                % volatile component base layer [wt H2O]
dv       =  1e-6;                % amplitude of random noise [wt H2O]

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
fin      =  1;                   % ingassing factor (0 = no ingassing; 1 = free flow ingassing)
fout     =  1;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)
tau_T    =  8*hr;                % wall cooling/assimilation time [s]
tau_a    =  4*hr;                % wall cooling/assimilation tie [s]
Twall    =  300;                 % wall temperature [degC] (nan = insulating)
cwall    =  nan;                 % wall major component [wt SiO2] (nan = no assimilation)
vwall    =  nan;                 % wall volatile component [wt H2O] (nan = no assimilation)
tewall   =  [nan,nan,nan,nan];   % wall trace elements [wt ppm] (nan = no assimilation)
irwall   =  [nan,nan,nan,nan];   % wall isotope ratios [delta] (nan = no assimilation)

% set thermo-chemical material parameters
calID    =  'morb';              % phase diagram calibration
kT0      =  5;                   % thermal conductivity [W/m/K]
cP       =  1200;                % heat capacity [J/kg/K]
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  400;                 % entropy change of exsolution [J/kg]

% set model buoyancy parameters
d0       =  1e-3;                % crystal/bubble size [m]
g0       =  10.;                 % gravity [m/s2]

% set numerical model parameters
CFL      =  0.90;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
BCA      =  {'',''};             % boundary condition on advection (top/bot, sides)
rtol     =  1e-3;                % outer its relative tolerance
atol     =  1e-6;                % outer its absolute tolerance
maxit    =  10;                  % maximum outer its
lambda   =  0.25;                % iterative lag parameter equilibration
etareg   =  1e0;                 % viscosity regularisation parameter


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

