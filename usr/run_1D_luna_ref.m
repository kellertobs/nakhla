% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '1D_luna_ref';       % run identifier
opdir    =  '../out';            % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  500;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  1000e3;              % chamber depth [m]
N        =  2300;                 % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  h;                   % chamber width [m]

% set model timing parameters
Nt       =  5e5;                 % number of time steps to take
tend     =  100*yr;              % end time for simulation [s]
dt       =  1*hr;                % initial time step [s]
dtmax    =  1*yr;                % maximum time step [s]

% set initial thermo-chemical state
T0       =  1705;                % temperature top  layer [deg C]
T1       =  T0;                  % temperature base layer [deg C]
c0       =  [0.36,0.31,0.32,0.01,0.0];   % components (maj comp, H2O) top layer [wt] (will be normalised to unit sum!)
c1       =  c0;                          % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)
dcr      =  [1/2,1/2,-1/2,-1/2,0]*0e-5;  % amplitude of random noise [wt]
zlay     =  2.0;                 % layer thickness (relative to domain depth D)

% set thermo-chemical boundary parameters
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  h;                   % boundary layer width [m]
tau_T    =  1*yr;                % wall cooling/assimilation time [s]
Twall    =  [0,1900,nan];        % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
Ptop     =  1e5;                 % top pressure [Pa]

% set thermo-chemical material parameters
calID    =  'luna';              % phase diagram calibration
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  400;                 % entropy change of exsolution [J/kg]
aT       =  5e-5;                % thermal expansivity [1/K]
cP       =  1100;                % heat capacity [J/kg/K]

% set buoyancy parameters
g0       =  1.62;                % gravity [m/s2]
dx       =  1e-3;                % crystal size [m]
bPx      =  1e-11;               % solid compressibility [1/Pa]
bPm      =  3e-11;               % melt compressibility [1/Pa]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-9;                % outer its absolute tolerance
maxit    =  20;                  % maximum outer its
Delta    =  2*D/50;              % correlation length for eddy diffusivity
Prt      =  1;                   % turbulent Prandtl number (ratio of momentum to heat diffusivity)
mink     =  1e3;                 % minimum eddy diffusivity constant

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

