% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID     =  '2D_DEMO_rchrg_inv';     % run identifier
restart   =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop       =  10;                  % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  0;                   % switch on to save output to file
plot_cv   =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D         =  20;                  % chamber depth [m]
N         =  100;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  D/1;                 % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt        =  5e5;                 % number of time steps to take
tend      =  1*yr;                % end time for simulation [s]
dt        =  1;                  % initial time step [s]

% set initial thermo-chemical state
zlay      =  0.5;                 % layer thickness (relative to domain depth D)
dlay      =  0.05;                 % random perturbation to layer thickness (relative to grid spacing h)
wlay_T    =  h/D/1;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c    =  h/D/1;               % thickness of smooth layer boundary (relative to domain depth D)
smth      =  15;                  % regularisation of initial random perturbation
T1        =  675;                 % temperature top  layer [deg C]
T0        =  1210;                % temperature base layer [deg C]
c1        =  [0.1  0.1  0.3  6.5  38  55  3.0]/100; % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)
c0        =  [11   17   35   31    3   3  0.5]/100; % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
dcr       =  [1,1,1,-1,-1,-1,0]*1e-6;
trc1      =  [0.1,0.3,0.5,3,10,2];  % trace elements top layer [wt ppm]
trc0      =  [1,1,1,1,1,1];         % trace elements base layer [wt ppm]
dr_trc    =  [1,1,1,-1,-1,-1]*1e-6; % trace elements random noise

% set thermo-chemical boundary parameters
periodic  =  1;
bndmode   =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w     =  h;                   % boundary layer width [m]
tau_T     =  bnd_w^2/1e-6;        % wall cooling/assimilation time [s]
Twall     =  [nan,T1,nan];        % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall     =  nan(3,6);
Ptop      =  1.5e8;               % top pressure [Pa]
fin       =  0;
fout      =  1;

% set thermo-chemical material parameters
calID     =  'DEMO';              % phase diagram calibration

% set numerical model parameters
TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL       =  0.75;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-4;                % outer its relative tolerance
atol      =  1e-8;                % outer its absolute tolerance
maxit     =  20;                  % maximum outer its
Delta_cnv =  h;                   % correlation length for eddy, convection diffusivity (multiple of h, 0.5-1)
Delta_sgr =  dx0*10;              % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)
alpha     =  0.75;
etamin    =  100;

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

