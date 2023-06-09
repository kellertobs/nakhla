% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '2D_rchrg';          % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  100;                 % chamber depth [m]
N        =  150;                 % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  D/2;                 % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt       =  5e5;                 % number of time steps to take
tend     =  1*yr;                % end time for simulation [s]
dt       =  1;

% set initial thermo-chemical state
T0       =  715;                 % temperature top layer [deg C]
T1       =  1200;                % temperature top layer [deg C]
c0       =  [0.01,0.01,0.28,0.70,0.04];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1       =  [0.14,0.49,0.28,0.09,0.02];  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr      =  [1/2,1/2,-1/3,-1/3,-1/3]*1e-4;
zlay     =  0.9;                 % layer thickness (relative to domain depth D)
dlay     =  0.0;                 % random perturbation to layer thickness (relative to grid spacing h)
wlay_T   =  1.0*h/D;             % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  0.5*h/D;             % thickness of smooth layer boundary (relative to domain depth D)

% set model trace and isotope geochemistry parameters (must match # trace elements and isotope ratios in calibration!)
te0      =  [1 ,1,  1,   1];     % trace elements top layer [wt ppm]
te1      =  [10,3,0.1,0.01];     % trace elements top layer [wt ppm]
ir0      =  [-1,5];              % isotope ratios top layer [delta]
ir1      =  [ 1,1];              % isotope ratios top layer [delta]

% set thermo-chemical boundary parameters
Ptop     =  1.25e8;              % top pressure [Pa]
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  h;                   % boundary layer width [m]
bnd_h    =  [0,0,0];             % internal wall rock layer thickness [m]
fin      =  0;                   % ingassing factor (0 = no ingassing; 1 = free flow ingassing)
fout     =  1;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)
tau_T    =  8*hr;                % wall cooling/assimilation time [s]
tau_a    =  8*hr;                % wall cooling/assimilation tie [s]
Twall    =  [300,1200,nan];      % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  [nan,nan,nan,nan,nan; ...
             0.14,0.49,0.28,0.09,0.02; ...
             nan,nan,nan,nan,nan]; % [top,bot,sds] wall rock major component [wt SiO2] (nan = no assimilation)
tewall   =  [nan,nan,nan,nan; ...
              10,  3,0.1,0.01; ...
             nan,nan,nan,nan];   % [top,bot,sds] wall rock trace elements [wt ppm] (nan = no assimilation)
irwall   =  [nan,nan; ...
               1,  1; ...
             nan,nan];           % [top,bot,sds] wall rock isotope ratios [delta] (nan = no assimilation)

% set thermo-chemical material parameters
calID    =  'andesSVZ';          % phase diagram calibration
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  400;                 % entropy change of exsolution [J/kg]

% set model buoyancy parameters
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.75;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-9;                % outer its absolute tolerance
maxit    =  20;                  % maximum outer its
Delta    =  2*D/50;              % correlation length for eddy viscosity


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

