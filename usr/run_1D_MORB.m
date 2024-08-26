% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '1D_MORB';           % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  10;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  0;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  10;                  % chamber depth [m]
N        =  100;                 % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  h;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt       =  5e5;                 % number of time steps to take
tend     =  1*yr;                % end time for simulation [s]
dt       =  36;                  % initial time step [s]

% set initial thermo-chemical state
T0       =  1175;                % temperature top  layer [deg C]
T1       =  T0;                  % temperature base layer [deg C]
c0       =  [0.03  0.22  0.32  0.31  0.09  0.03  0.005];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1       =  c0;                  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr      =  [0,0,0,0,0,0,0];
dcg      =  [0,0,0,0,0,0,0];

% set thermo-chemical boundary parameters
periodic =  1;
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  0.1;                 % boundary layer width [m]
tau_T    =  (2*bnd_w)^2/1e-6;    % wall cooling/assimilation time [s]
Twall    =  [300,300,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  nan(3,7);
Ptop     =  2.0e8;               % top pressure [Pa]
fin      =  0;
fout     =  1;

% set thermo-chemical material parameters
calID    =  'MORB';              % phase diagram calibration
% aTm      =  4e-5;                % melt  thermal expansivity [1/K]
% aTx      =  4e-5;                % xtal  thermal expansivity [1/K]
% aTf      =  4e-4;                % fluid thermal expansivity [1/K]
% kTm      =  4;                   % melt  thermal conductivity [W/m/K]
% kTx      =  4;                   % xtal  thermal conductivity [W/m/K]
% kTf      =  4;                   % fluid thermal conductivity [W/m/K]
% cPm      =  1200;                % melt  heat capacity [J/kg/K]
% cPx      =  1200;                % xtal  heat capacity [J/kg/K]
% cPf      =  1200;                % fluid heat capacity [J/kg/K]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-8;                % outer its absolute tolerance
maxit    =  50;                  % maximum outer its
Delta_cnv=  h/2;                 % correlation length for eddy, convection diffusivity (multiple of h, 0.5-1)
Delta_sgr=  dx0*10;              % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

