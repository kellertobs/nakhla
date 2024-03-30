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
N        =  120;                 % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  D/1;                 % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt       =  5e5;                 % number of time steps to take
tend     =  1*yr;                % end time for simulation [s]
dt       =  24;                  % initial time step [s]

% set initial thermo-chemical state
zlay     =  0.80;                 % layer thickness (relative to domain depth D)
wlay_T   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0       =  720;                 % temperature top  layer [deg C]
T1       =  1245;                % temperature base layer [deg C]
c0       =  [0.01  0.02  0.04  0.08  0.16  0.69  0.03];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1       =  [0.09  0.12  0.30  0.35  0.06  0.08  0.003]; % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr      =  [1,1,1,-1,-1,-1,0]*1e-4;
trc0     =  [0.1,0.3,0.5,3,10,2];       % trace elements top layer [wt ppm]
trc1     =  [1,1,1,1,1,1];       % trace elements base layer [wt ppm]
dr_trc   =  [0,0,0,0,0,0];       % trace elements random noise

% set thermo-chemical boundary parameters
periodic =  0;
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  h;                   % boundary layer width [m]
tau_T    =  4*hr;                % wall cooling/assimilation time [s]
Twall    =  [nan,T0,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  nan(3,7);
Ptop     =  2.0e8;               % top pressure [Pa]
fin      =  0;
fout     =  1;

% set thermo-chemical material parameters
calID    =  'MORB';              % phase diagram calibration

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.75;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-3;                % outer its relative tolerance
atol     =  1e-8;                % outer its absolute tolerance
maxit    =  20;                  % maximum outer its
etamin   =  10;
gamma    =  1e-2;

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

