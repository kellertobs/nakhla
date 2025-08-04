 % prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID     =  '1D_DEMO';           % run identifier
restart   =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nrh       =  1e2;                 % record diagnostic history every 'nrh' time steps
nop       =  1e3;                 % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  1;                   % switch on to save output to file
colourmap = 'lapaz';              % choose colourmap ('ocean','lipari','lajolla','lapaz','navia','batlow(W/K)','glasgow')

% set model domain parameters
D         =  10;                  % chamber depth [m]
N         =  200;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  h;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt        =  5e5;                 % number of time steps to take
tend      =  1*yr;                % end time for simulation [s]
dt        =  1;                   % initial time step [s]

% set initial thermo-chemical state
init_mode =  'liquidus';
T0        =  -5;                  % initial temperature [deg C]
c0        =  [10  19  32  30  5  4  0.5]/100;  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
dcr       =  [1,1,1,-1,-1,-1,0]*1e-4;
dr_trc    =  [1,1,1,-1,-1,-1  ]*1e-4; % trace elements random noise

% set thermo-chemical boundary parameters
periodic  =  1;                   % periodic side boundaries
bndmode   =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w     =  h;                   % boundary layer width [m]
tau_T     =  1*hr;                % wall cooling/assimilation time [s]
Twall     =  [300,300,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
Ptop      =  1.5e8;               % top pressure [Pa]

% set thermo-chemical material parameters
calID     =  'DEMO';              % phase diagram calibration

% set numerical model parameters
TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL       =  1.00;                % (physical) time stepping courant number (multiplies stable step) [0,1]
alpha     =  0.75;                % iterative step size parameter
rtol      =  1e-4;                % outer its relative tolerance
atol      =  1e-9;                % outer its absolute tolerance
maxit     =  20;                  % maximum outer its
Delta_cnv =  D/10;                % correlation length for eddy, convection diffusivity (multiple of h, 0.5-1)
Delta_sgr =  dx0*10;              % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)
kmax      =  1e-2;                % maximum diffusivity
Prt       =  3;                   % turbulent Prandtl number (ratio of momentum to heat diffusivity)
Sct       =  3;                   % turbulent Schmidt number (ratio of momentum to mass diffusivity)


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

