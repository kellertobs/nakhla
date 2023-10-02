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
N        =  100;                 % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  D;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt       =  5e5;                 % number of time steps to take
tend     =  1*yr;                % end time for simulation [s]
dt       =  10;

% set initial thermo-chemical state
zlay     =  0.9;                 % layer thickness (relative to domain depth D)
wlay_T   =  1/D;                 % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  1/D;                 % thickness of smooth layer boundary (relative to domain depth D)
T0       =  735;                 % temperature top  layer [deg C]
T1       =  1215;                % temperature base layer [deg C]
c0       =  [0.01  0.01  0.03  0.25  0.70  0.050];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1       =  [0.13  0.20  0.53  0.11  0.03  0.005];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
dcr      =  [-1,-1,-1,1,1,1]*1e-6;
dcg      =  [0,0,0,0,0,0];
trc0     =  [0.1,0.3,2,3  ,10 ,1];  % trace elements top layer [wt ppm]
trc1     =  [10 ,3  ,1,0.3,0.1,2];  % trace elements base layer [wt ppm]

% set thermo-chemical boundary parameters
periodic =  1;
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  h;                   % boundary layer width [m]
tau_T    =  24*hr;               % wall cooling/assimilation time [s]
Twall    =  [T0,T0,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  [nan,nan,nan,nan,nan,nan; ...
             nan,nan,nan,nan,nan,nan; ...
             nan,nan,nan,nan,nan,nan];
Ptop     =  1.5e8;               % top pressure [Pa]
fin      =  0;
fout     =  1;

% set thermo-chemical material parameters
calID    =  'MORB';              % phase diagram calibration

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-3;                % outer its relative tolerance
atol     =  1e-8;                % outer its absolute tolerance
maxit    =  30;                  % maximum outer its
Delta    =  2*D/100;             % correlation length for eddy viscosity
etamin   =  100;

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

