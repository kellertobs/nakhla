% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '2D_fract_anh';      % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  120 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  5e5;                 % number of time steps to take
tend     =  1*yr;                % end time for simulation [s]

% set initial thermo-chemical state
T0       =  1250;                % temperature top layer [deg C]
c0       =  0.51;                % major component top layer [wt SiO2]
dcr      =  1e-4;                % amplitude of random noise [wt SiO2]
v0       =  0.00;                % volatile component top layer [wt H2O]
dvr      =  0e-4;                % amplitude of random noise [wt H2O]

% set model trace and isotope geochemistry parameters (must match # trace elements and isotope ratios in calibration!)
te0      =  [1,1,1,1];           % trace elements top layer [wt ppm]
ir0      =  [0, 1];              % isotope ratios top layer [delta]
dirr     =  [1, 0];              % isotope ratios random noise [delta]

% set thermo-chemical boundary parameters
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  h;                   % boundary layer width [m]
tau_T    =  6*hr;                % wall cooling/assimilation time [s]
Twall    =  [300,300,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
Ptop     =  1.25e8;              % top pressure [Pa]

% set thermo-chemical material parameters
calID    =  'default';           % phase diagram calibration
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  400;                 % entropy change of exsolution [J/kg]
tau_r    =  0;                   % reaction time scale (set to zero for quasi-equilibrium mode)

% set model buoyancy parameters
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-8;                % outer its absolute tolerance
maxit    =  20;                  % maximum outer its
cnvreg   =  10;                  % convection regularisation parameter
dtmax    =  10;                  % maximum time step [s]


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

