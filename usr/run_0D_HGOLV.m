% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '0D_HGOLV';     % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  1;                   % chamber depth [m]
L        =  D;                   % chamber width [m]
N        =  1;                   % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  1e4;                 % number of time steps to take
tend     =  1e4*hr;              % end time for simulation [s]
dt       =  hr/3;                % initial time step [s]
dtmax    =  hr/3;                % maximum time step [s]
 
% set initial thermo-chemical state
T0       =  1800;                % temperature top  layer [deg C]
T1       =  T0;                  % temperature base layer [deg C]
c0       =  [30  23  16  18  13  0]/100;  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1       =  c0;                  % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)
dcr      =  [0,0,0,0,0,0,0,0];
dcg      =  [0,0,0,0,0,0,0,0];

% set thermo-chemical boundary parameters
fractxtl =  1;                   % fractional crystallisation mode for 0-D (Nz=Nx=1)
fractmlt =  0;                   % fractional melting mode for 0-D (Nz=Nx=1)
fractres =  0.01;                % residual fraction for fractionation mode
dPdT     =  4.75e6;                 % decompression rate for 0D models
bndmode  =  1;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  1e16;                % boundary layer width [m]
tau_T    =  D^2/1e-6;            % wall cooling/assimilation time [s]
Twall    =  [0,0,nan];           % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
Ptop     =  5e9;                 % top pressure [Pa]
fin      =  0;                   % ingassing factor (0 = no ingassing; 1 = free flow ingassing)
fout     =  0;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)

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
dx0      =  1e-2;                % crystal size [m]
dm0      =  1e-2;                % melt film size [m]
bPx      =  1e-11;               % solid compressibility [1/Pa]
bPm      =  4e-11;               % melt compressibility [1/Pa]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  1.00;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-5;                % outer its relative tolerance
atol     =  1e-8;                % outer its absolute tolerance
maxit    =  50;                  % maximum outer its
Pcouple  =  1;                   % coupling phase equilibria and material properties to dynamic pressure
alpha    =  0.75;

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

