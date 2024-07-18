% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '0D_TORF_fractxtl25';% run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  0.1;                 % chamber depth [m]
L        =  0.1;                 % chamber width [m]
N        =  1;                   % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  1e4;                 % number of time steps to take
tend     =  12*hr;               % end time for simulation [s]
dt       =  12;                  % initial time step [s]
dtmax    =  12;                  % maximum time step [s]

% set initial thermo-chemical state
T0       =  1300;                % temperature top  layer [deg C]
T1       =  T0;                  % temperature base layer [deg C]
c0       =  [0.034  0.413  0.412  0.130  0.010  0.005];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1       =  c0;                  % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)
dcr      =  [0,0,0,0,0,0,0,0];
dcg      =  [0,0,0,0,0,0,0,0];

% set thermo-chemical boundary parameters
fractxtl =  1;                   % fractional crystallisation mode for 0-D (Nz=Nx=1)
fractmlt =  0;                   % fractional melting mode for 0-D (Nz=Nx=1)
fractres =  0.25;                % residual fraction for fractionation mode
bndmode  =  1;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  1e16;                % boundary layer width [m]
tau_T    =  D^2/1e-6;            % wall cooling/assimilation time [s]
Twall    =  [300,300,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  nan(3,6);
Ptop     =  1.5e8;               % top pressure [Pa]
fin      =  0;                   % ingassing factor (0 = no ingassing; 1 = free flow ingassing)
fout     =  0;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)

% set thermo-chemical material parameters
calID    =  'TORF';              % phase diagram calibration

% set chamber pressure parameters
Pchmb0   =  0;                  % initial chamber pressure [Pa]
eta_wall =  1;                  % wall rock viscosity [Pas]
mod_wall =  0;                  % wall rock elastic modulus [Pa]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  1.00;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-9;                % outer its absolute tolerance
maxit    =  50;                  % maximum outer its


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

