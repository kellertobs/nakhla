% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '0D_fract_sat';      % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  1;                   % chamber depth [m]
L        =  1;                   % chamber width [m]
N        =  1 + 2;               % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  2000;                % number of time steps to take
tend     =  14*hr;               % end time for simulation [s]
dt       =  36;                  % initial time step [s]
dtmax    =  36;                  % maximum time step [s]

% set initial thermo-chemical state
T0       =  1220;                % temperature top layer [deg C]
c0       =  0.51;                % major component top layer [wt SiO2]
v0       =  0.04;                % volatile component top layer [wt H2O]

% set thermo-chemical boundary parameters
Ptop     =  1.25e8;              % top pressure [Pa]
bnd_w    =  D;                   % boundary layer width [m]
tau_T    =  12*hr;               % wall cooling/assimilation time [s]
Twall    =  300;                 % wall temperature [degC] (nan = insulating)

% set thermo-chemical material parameters
calID    =  'default';           % phase diagram calibration
kT0      =  5;                   % thermal conductivity [W/m/K]
cP       =  1200;                % heat capacity [J/kg/K]
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  400;                 % entropy change of exsolution [J/kg]

% set numerical model parameters
CFL      =  0.75;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-9;                % outer its absolute tolerance
maxit    =  50;                  % maximum outer its
lambda   =  0.50;                % iterative lag parameter equilibration


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

