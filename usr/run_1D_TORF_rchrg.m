% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID     =  '1D_TORF_rchrg_bot'; % run identifier
restart   =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop       =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  1;                   % switch on to save output to file
plot_cv   =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D         =  50;                  % chamber depth [m]
N         =  300;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  h;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt        =  5e5;                 % number of time steps to take
tend      =  1*yr;                % end time for simulation [s]
dt        =  48;                  % initial time step [s]

% set initial thermo-chemical state
zlay      =  0.9;                 % layer thickness (relative to domain depth D)
wlay_T    =  0*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c    =  0*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0        =  700;                 % temperature top  layer [deg C]
T1        =  1200;                % temperature base layer [deg C]
c0        =  [0.001 0.001 0.020 0.277 0.700 0.023];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1        =  [0.034 0.413 0.412 0.130 0.010 0.005];  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr       =  [1,1,1,-1,-1,-1,0]*0e-6;
trc0      =  [0.1,0.3,0.5,3,10,2];       % trace elements top layer [wt ppm]
trc1      =  [1,1,1,1,1,1];       % trace elements base layer [wt ppm]
dr_trc    =  [0,0,0,0,0,0];       % trace elements random noise

% set thermo-chemical boundary parameters
periodic  =  1;
bndmode   =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w     =  h;                   % boundary layer width [m]
tau_T     =  bnd_w^2/1e-6;        % wall cooling/assimilation time [s]
Twall     =  [nan,T0,nan];        % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall     =  nan(3,6);
Ptop      =  1.5e8;               % top pressure [Pa]
fin       =  0;
fout      =  1;

% set thermo-chemical material parameters
calID     =  'TORF';              % phase diagram calibration

% set numerical model parameters
TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL       =  0.25;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-3;                % outer its relative tolerance
atol      =  1e-8;                % outer its absolute tolerance
maxit     =  15;                  % maximum outer its
Delta_cnv =  h/2;                 % correlation length for eddy, convection diffusivity (multiple of h, 0.5-1)
Delta_sgr =  h/20;                % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)
Sct       =  6;
Prt       =  3;
alpha     =  0.5;
beta      =  0;

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

