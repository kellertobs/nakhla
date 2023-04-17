% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '0D_andesSVZ_hyd';   % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  50;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  0;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  1;                   % chamber depth [m]
L        =  1;                   % chamber width [m]
N        =  1;                   % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  900;                 % number of time steps to take
tend     =  9*hr;                % end time for simulation [s]
dt       =  36;                  % initial time step [s]
dtmax    =  36;                  % maximum time step [s]

% set initial thermo-chemical state
T0       =  1190;                % temperature top layer [deg C]
c0       =  [0.10,0.28,0.59,0.03,0.02]; % components (maj comp, H2O) top layer [wt] (will be normalised to unit sum!)
c1       =  c0;                         % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)

% set thermo-chemical boundary parameters
bndmode  =  1;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  1e16;                % boundary layer width [m]
tau_T    =  8*hr;                % wall cooling/assimilation time [s]
Twall    =  [300,300,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  [nan,nan,nan,nan,nan; ...
             nan,nan,nan,nan,nan; ...
             nan,nan,nan,nan,nan];
Ptop     =  1.25e8;              % top pressure [Pa]
fin      =  0;                   % ingassing factor (0 = no ingassing; 1 = free flow ingassing)
fout     =  0;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)

% set thermo-chemical material parameters
calID    =  'andesSVZ';          % phase diagram calibration
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  400;                 % entropy change of exsolution [J/kg]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  1.00;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-5;                % outer its relative tolerance
atol     =  1e-9;                % outer its absolute tolerance
maxit    =  50;                  % maximum outer its


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

