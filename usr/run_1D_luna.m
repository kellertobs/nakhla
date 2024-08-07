% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '1D_luna';           % run identifier
opdir    =  '../out';            % output directory
restart  = -1;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  200;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  1000e3;              % chamber depth [m]
N        =  300;                 % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  h;                 % chamber width [m]

% set model timing parameters
Nt       =  1e6;                 % number of time steps to take
tend     =  1e4*yr;              % end time for simulation [s]
dt       =  1*hr;                % initial time step [s]
dtmax    =  1*yr;                % maximum time step [s]

% set initial thermo-chemical state
Tinit    = 'linear';             % T initial condition mode ('layer' or 'linear')
T0       =  1725;                % temperature top  layer [deg C]
T1       =  1725;                % temperature base layer [deg C]
c0       =  [0.30  0.31  0.10  0.20  0.05  0.04  0.00];  % components (maj comp, H2O) top layer [wt] (will be normalised to unit sum!)
c1       =  c0;                             % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)
dcr      =  [1,1,1,-1,-1,-1,0]*1e-4;  % amplitude of random noise [wt]
dcg      =  [0,0,0,0,0,0,0];          % amplitude of gaussian perturbation [wt]
zlay     =  2.0;                 % layer thickness (relative to domain depth D)

% set thermo-chemical boundary parameters
periodic =  1;
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  h;                   % boundary layer width [m]
tau_T    =  2*yr;                % wall cooling/assimilation time [s]
Twall    =  [0,nan,nan];         % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  nan(3,7);
Ptop     =  1e5;                 % top pressure [Pa]

% set thermo-chemical material parameters
calID    =  'LUNA';              % phase diagram calibration
aT       =  3e-5;                % thermal expansivity [1/K]
cP       =  1100;                % heat capacity [J/kg/K]

% set buoyancy parameters
g0       =  1.62;                % gravity [m/s2]
dx0      =  3e-3;                % crystal size [m]
dm0      =  3e-3;                % melt film size [m]
bPx      =  1e-11;               % solid compressibility [1/Pa]
bPm      =  2e-11;               % melt compressibility [1/Pa]

% switch off chamber pressure parameters
Pchmb0   =  0;                   % initial chamber pressure [Pa]
eta_wall =  1;                   % wall rock viscosity [Pas]
mod_wall =  0;                   % wall rock elastic modulus [Pa]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-8;                % outer its absolute tolerance
maxit    =  20;                  % maximum outer its
Delta_cnv=  h/2;                 % correlation length for eddy diffusivity
Delta_sgr=  h/2;                 % correlation length for phase fluctuation diffusivity

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

