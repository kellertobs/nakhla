% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID     =  '2D_TORF';           % run identifier
restart   = -1;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop       =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  1;                   % switch on to save output to file
plot_cv   =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D         =  10;                  % chamber depth [m]
N         =  150;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  D/2;                 % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt        =  5e5;                 % number of time steps to take
tend      =  1*yr;                % end time for simulation [s]
dt        =  36;                  % initial time step [s]

% set initial thermo-chemical state
T0        =  1200;                % temperature top  layer [deg C]
T1        =  T0;                  % temperature base layer [deg C]
c0        =  [0.034  0.413  0.412  0.130  0.010  0.005];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1        =  c0;                  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr       =  [1,1,1,-1,-1,-1,0]*1e-4;
dr_trc    =  [0,0,1,0,0,-1];      % trace elements random noise

% set thermo-chemical boundary parameters
periodic  =  1;
bndmode   =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w     =  0.1;                 % boundary layer width [m]
tau_T     =  bnd_w^2/1e-6;        % wall cooling/assimilation time [s]
Twall     =  [300,300,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall     =  nan(3,6);
Ptop      =  1.5e8;               % top pressure [Pa]
fin       =  0;
fout      =  1;

% set thermo-chemical material parameters
calID     =  'TORF';              % phase diagram calibration

% set segregation parameters
dm0       =  1e-3;                % melt film size [m]
dx0       =  1e-3;                % crystal size [m]
df0       =  1e-3;                % bubble size [m]

% set numerical model parameters
TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL       =  0.75;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-3;                % outer its relative tolerance
atol      =  1e-8;                % outer its absolute tolerance
maxit     =  15;                  % maximum outer its
Delta_cnv =  h/2;                 % correlation length for eddy, convection diffusivity (multiple of h, 0.5-1)
Delta_sgr =  h/10;                % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

