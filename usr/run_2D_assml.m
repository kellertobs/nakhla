% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '2D_assml';          % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  10;                  % chamber depth [m]
N        =  100;                 % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  D;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt       =  5e5;                 % number of time steps to take
tend     =  1*yr;                % end time for simulation [s]

% set initial thermo-chemical state
T0       =  1175;                % temperature top layer [deg C]
c0       =  0.51;                % major component top layer [wt SiO2]
dcr      =  1e-4;                % amplitude of random noise [wt SiO2]
v0       =  0.02;                % volatile component top layer [wt H2O]
dvr      =  1e-4;                % amplitude of random noise [wt H2O]

% set model trace and isotope geochemistry parameters (must match # trace elements and isotope ratios in calibration!)
te0      =  [1,1,1,1];           % trace elements top layer [wt ppm]
ir0      =  [-1,5];              % isotope ratios top layer [delta]

% set thermo-chemical boundary parameters
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_h    =  [10*h,10*h,0];       % boundary layer width [m]
bnd_w    =  1.5*h;               % boundary layer width [m]
tau_T    =  12*hr;               % wall cooling/assimilation time [s]
tau_a    =  24*hr;               % wall cooling/assimilation time [s]
Twall    =  [500,500,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  [0.70,0.70,nan];     % [top,bot,sds] wall rock major component [wt SiO2] (nan = no assimilation)
vwall    =  [0.04,0.04,nan];     % [top,bot,sds] wall rock volatile component [wt H2O] (nan = no assimilation)
tewall   =  [0.01,0.1,3,10; ...
             0.01,0.1,3,10; ...
             nan,nan,nan,nan];   % [top,bot,sds] wall rock trace elements [wt ppm] (nan = no assimilation)
irwall   =  [1,1; ...
             1,1; ...
             nan,nan];           % [top,bot,sds] wall rock isotope ratios [delta] (nan = no assimilation)
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
