% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '1D_luna_newpar';      % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
save_op  =  0;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  1350e3;              % chamber depth [m]
N        =  500;                 % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  h;                   % chamber width [m]

% set model timing parameters
Nt       =  1e6;                 % number of time steps to take
tend     =  1e4*yr;              % end time for simulation [s]
dt       =  1*hr;                % initial time step [s]
dtmax    =  1*yr;                % maximum time step [s]

% set initial thermo-chemical state
init_mode= 'layer';            % T initial condition mode ('layer' or 'linear')
T0       =  1750;                % temperature top  layer [deg C]
T1       =  1750;                % temperature base layer [deg C]
c0       =  [0.30  0.31  0.10  0.20  0.05  0.04  0.00];  % components (maj comp, H2O) top layer [wt] (will be normalised to unit sum!)
c1       =  [0.30  0.31  0.10  0.20  0.05  0.04  0.00];                             % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)
dcr      =  [1,1,1,-1,-1,-1,0]*0e-4;  % amplitude of random noise [wt]
dcg      =  [0,0,0,0,0,0,0];          % amplitude of gaussian perturbation [wt]
zlay     =  2.0;                 % layer thickness (relative to domain depth D)

% set thermo-chemical boundary parameters
periodic =  1;
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  D/200;               % boundary layer width [m]
tau_T    =  yr/10;               % wall cooling/assimilation time [s]
Twall    =  [0,nan,nan];         % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  nan(3,7);
Ptop     =  1e5;                 % top pressure [Pa]

% set thermo-chemical material parameters
calID    =  'LUNA';              % phase diagram calibration
aTm      =  4.5e-5;              % melt  thermal expansivity [1/K]
aTx      =  2e-5;                % xtal  thermal expansivity [1/K]
aTf      =  1e-4;                % fluid thermal expansivity [1/K]
kTm      =  1;                   % melt  thermal conductivity [W/m/K]
kTx      =  4;                   % xtal  thermal conductivity [W/m/K]
kTf      =  0.5;                 % fluid thermal conductivity [W/m/K]
cPm      =  1300;                % melt  heat capacity [J/kg/K]
cPx      =  1000;                % xtal  heat capacity [J/kg/K]
cPf      =  2000;                % fluid heat capacity [J/kg/K]

% set buoyancy parameters
g0       =  1.62;                % gravity [m/s2]
dx0      =  1e-2;                % crystal size [m]
dm0      =  1e-2;                % melt film size [m]
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
atol     =  1e-7;                % outer its absolute tolerance
maxit    =  15;                  % maximum outer its
gamma    =  0.001;               % horizontal drag
Delta_cnv=  h;                   % correlation length for eddy, convection diffusivity (multiple of h, 0.5-1)
Delta_sgr=  dx0*10;              % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)
r_unst   =  0.0;                 % coefficient boosting convecive mixing diffusivity where unstable density gradient present
alpha    =  0.5;
beta     =  0.0;

%*****  RUN NAKHLA MODEL  *************************************************
run([srcdir,'/main'])
%**************************************************************************

