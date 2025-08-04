% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '2D_HGOLV_slowcool';     % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  10;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  0;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  450e3;               % chamber depth [m]
N        =  100;                 % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  D/1;                 % chamber width [m]

% set model timing parameters
Nt       =  1e6;                 % number of time steps to take
tend     =  1e4*yr;              % end time for simulation [s]
dt       =  0.1*hr;                % initial time step [s]
dtmax    =  10*hr;

% set initial thermo-chemical state
init_mode= 'linear';
T0       =  1780;                % temperature top  layer [deg C]
T1       =  1815;                  % temperature base layer [deg C]
c0       =  [30  23  16  18  13  0]/100;  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1       =  c0;                  % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)
dcr      =  [1,1,-2/3,-2/3,-2/3,0]*1e-3;
dcg      =  [0,0,0,0,0,0,0,0];
dTr      =  2;

% set thermo-chemical boundary parameters
periodic =  1;
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  h;                   % boundary layer width [m]
tau_T    =  10*yr/3;             % wall cooling/assimilation time [s]
Twall    =  [0,1950,nan];        % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
Ptop     =  1e5;                 % top pressure [Pa]

% set thermo-chemical material parameters
calID    =  'HGOLV';              % phase diagram calibration
aTm      =  4e-5;                % melt  thermal expansivity [1/K]
aTx      =  1e-5;                % xtal  thermal expansivity [1/K]
kTm      =  1;                   % melt  thermal conductivity [W/m/K]
kTx      =  4;                   % xtal  thermal conductivity [W/m/K]
cPm      =  1200;                % melt  heat capacity [J/kg/K]
cPx      =  1000;                % xtal  heat capacity [J/kg/K]

% set buoyancy parameters
g0       =  3.7;
dx0      =  1e-3;                % crystal size [m]
dm0      =  1e-3;                % melt film size [m]
bPx      =  1e-11;               % solid compressibility [1/Pa]
bPm      =  4e-11;               % melt compressibility [1/Pa]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.75;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-8;                % outer its absolute tolerance
maxit    =  20;                  % maximum outer its
gamma    =  0.01;                % horizontal drag
Delta_cnv=  h;                   % correlation length for eddy, convection diffusivity (multiple of h, 0.5-1)
Delta_sgr=  dx0*10;              % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)
alpha    =  0.75;                % iterative step size
Pcouple  =  0;
Prt      =  3;
Sct      =  3;
etamin   =  1e+4;                 % minimum viscosity
kmin     =  1e-9;                 % minimum eddy diffusivity
kmax     =  1e+6;                 % maximum eddy diffusivity
maxcmp   =  0.01;
lambda1  =  0e-13;

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

