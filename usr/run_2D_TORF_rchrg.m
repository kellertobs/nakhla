% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID     =  '2D_TORF_rchrg_bot_3'; % run identifier
srcdir    =  '../src';            % output directory
outdir    =  '../out';            % output directory
restart   =  0;                   % restart from file (0: new run; -1: restart from last; >1: restart from specified frame)
nop       =  10;                  % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  0;                   % switch on to save output to file

% set model domain parameters
D         =  20;                  % chamber depth [m]
N         =  100;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  D/2;                 % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt        =  1e6;                 % number of time steps to take
tend      =  10*yr;               % end time for simulation [s]
dt        =  10;                  % initial time step [s]

% set initial thermo-chemical state
zlay      =  0.50;                % layer thickness (relative to domain depth D)
dlay      =  0.01;                % random perturbation to layer thickness (relative to grid spacing h)
wlay_T    =  1*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c    =  1*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0        =  700;                 % temperature top  layer [deg C]
T1        =  1180;                % temperature base layer [deg C]
c0        =  [0.01 0.03 0.07 0.33 0.56 0.030];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1        =  [0.03 0.42 0.41 0.13 0.01 0.005];  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr       =  [1,1,1,-1,-1,-1,0]*1e-5;
trc0      =  [0.1,0.3,0.5,3,10,2];  % trace elements top layer [wt ppm]
trc1      =  [1,1,1,1,1,1];       % trace elements base layer [wt ppm]
dr_trc    =  [0,0,0,0,0,0];       % trace elements random noise

% set thermo-chemical boundary parameters
periodic  =  1;                   % set side boundaries to periodic
bndmode   =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w     =  h;                   % boundary layer width [m]
tau_T     =  bnd_w^2/1e-6;        % wall cooling/assimilation time [s]
Twall     =  [nan,T0,nan];        % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
Ptop      =  1.5e8;               % top pressure [Pa]

% set thermo-chemical material parameters
calID     =  'TORF';              % phase diagram calibration

% set numerical model parameters
TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL       =  0.75;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-3;                % outer its relative tolerance
atol      =  1e-8;                % outer its absolute tolerance
maxit     =  20;                  % maximum outer its
alpha     =  0.75;                % iterative step size parameter
beta      =  0.125;                % iterative damping parameter
gamma     =  1e-2;                % artificial horizontal inertia parameter (only applies if periodic)
Delta_cnv =  h/2;                 % correlation length for eddy diffusivity (multiple of h, 0.5-1)
Delta_sgr =  dx0*20;              % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)
Prt       =  3;                   % turbulent Prandtl number (ratio of momentum to heat diffusivity)
Sct       =  3;                   % turbulent Schmidt number (ratio of momentum to mass diffusivity)


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

