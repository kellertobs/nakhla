% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '2D_layer';          % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  0;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  100 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  5e5;                 % number of time steps to take
tend     =  1*yr;                % end time for simulation [s]

% set initial thermo-chemical state
T0       =  1035;                % temperature top layer [deg C]
T1       =  1150;                % temperature base layer [deg C]
c0       =  0.58;                % major component top layer [wt SiO2]
c1       =  0.51;                % major component base layer [wt SiO2]
dcr      =  1e-4;                % amplitude of random noise [wt SiO2]
v0       =  0.01;                % volatile component top layer [wt H2O]
v1       =  0.03;                % volatile component base layer [wt H2O]
dvr      =  0e-4;                % amplitude of random noise [wt H2O]
zlay     =  0.8;                 % layer thickness (relative to domain depth D)
dlay     =  0.25;                % random perturbation to layer thickness (relative to grid spacing h)
wlay_T   =  1.0*h/D;             % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  0.5*h/D;             % thickness of smooth layer boundary (relative to domain depth D)

% set model trace and isotope geochemistry parameters (must match # trace elements and isotope ratios in calibration!)
te0      =  [0.1,0.3,1,3];       % trace elements top layer [wt ppm]
te1      =  [3,1,0.3,0.1];       % trace elements base layer [wt ppm]
ir0      =  [-1, 1];             % isotope ratios top layer [delta]
ir1      =  [ 1,-1];             % isotope ratios base layer [delta]

% set thermo-chemical boundary parameters
bnd_w    =  h;                   % boundary layer width [m]
tau_T    =  12*hr;               % wall cooling/assimilation time [s]
Twall    =  300;                 % wall temperature [degC] (nan = insulating)
Ptop     =  1.25e8;              % top pressure [Pa]

% set thermo-chemical material parameters
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  400;                 % entropy change of exsolution [J/kg]

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
dtmax    =  5;                   % maximum time step [s]


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

