% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '2D_MORB_D200_np';        % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  50;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  200;                   % chamber depth [m]
N        =  90;                  % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  D*2/3;               % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt       =  5e5;                 % number of time steps to take
tend     =  1*yr;                % end time for simulation [s]
dt       =  36;                  % initial time step [s]

% set initial thermo-chemical state
T0       =  1200;                % temperature top  layer [deg C]
T1       =  T0;                  % temperature base layer [deg C]
c0       =  [0.08 0.24 0.41 0.24 0.03 0.005];  % components (maj comp, H2O) top layer [wt] (will be normalised to unit sum!)
c1       =  c0;                  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr      =  [-1,-1,-1,1,1,1]*1e-5;
dcg      =  [0,0,0,0,0,0];
ir0      =  [0, 0];              % isotope ratios random noise [delta]
dirr     =  [-1, 1];             % isotope ratios random noise [delta]

% set thermo-chemical boundary parameters
periodic =  0;
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  h;                   % boundary layer width [m]
tau_T    =  12*hr;               % wall cooling/assimilation time [s]
Twall    =  [300,300,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  [nan,nan,nan,nan,nan,nan; ...
             nan,nan,nan,nan,nan,nan; ...
             nan,nan,nan,nan,nan,nan];
Ptop     =  2.0e8;               % top pressure [Pa]
fin      =  0;
fout     =  1;

% set thermo-chemical material parameters
calID    =  'MORB';              % phase diagram calibration
Dsx      = -350;                 % entropy change of crystallisation [J/kg]
Dsf      =  450;                 % entropy change of exsolution [J/kg]
dx       =  0.001;

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.75;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-3;                % outer its relative tolerance
atol     =  1e-8;                % outer its absolute tolerance
maxit    =  30;                  % maximum outer its
Delta    =  2*h;                 % correlation length for eddy viscosity
etamin   =  0.01;

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

